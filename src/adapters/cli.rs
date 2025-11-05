use clap::Parser;
use std::fs;
use std::io::{self, Read};

use crate::error::AppError;
use crate::models::{Assumptions, Inputs};
use crate::salinity::calculator::CalculationSummary;

#[derive(Parser, Debug)]
#[command(author, version, about = "Salinity calculator (TEOS-10) — optional JSON output", long_about = None)]
pub struct Args {
    #[arg(long)]
    json: bool,
    #[arg(
        long,
        value_name = "FILE",
        help = "JSON file with inputs and optional assumptions; '-' reads from stdin"
    )]
    input: Option<String>,
    #[arg(
        long,
        value_name = "JSON",
        help = "Inline JSON for inputs (overrides --input)"
    )]
    inputs_json: Option<String>,
    #[arg(
        long,
        value_name = "JSON",
        help = "Inline JSON for assumptions (optional, supplements --inputs-json)"
    )]
    assumptions_json: Option<String>,
}

fn parse_inline_inputs(
    inputs_json: &str,
    assumptions_json: Option<&String>,
) -> Result<(Inputs, Assumptions), AppError> {
    let inputs: Inputs =
        serde_json::from_str(inputs_json).map_err(|source| AppError::ParseInputsJson { source })?;

    let assumptions = match assumptions_json {
        Some(s) => serde_json::from_str::<Assumptions>(s)
            .map_err(|source| AppError::ParseAssumptionsJson { source })?,
        None => Assumptions::default(),
    };

    Ok((inputs, assumptions))
}

fn parse_cmd_input_doc(doc: &str) -> Result<(Inputs, Assumptions), AppError> {
    let parsed: CmdInput =
        serde_json::from_str(doc).map_err(|source| AppError::ParseCmdInputJson { source })?;
    Ok((parsed.inputs, parsed.assumptions.unwrap_or_default()))
}

pub fn parse_inputs(args: &Args) -> Result<(Inputs, Assumptions), AppError> {
    match (&args.inputs_json, &args.input) {
        (Some(inputs_json), _) => parse_inline_inputs(inputs_json, args.assumptions_json.as_ref()),
        (None, Some(path)) if path == "-" => {
            let mut s = String::new();
            io::stdin()
                .read_to_string(&mut s)
                .map_err(|source| AppError::ReadStdin { source })?;
            parse_cmd_input_doc(&s)
        }
        (None, Some(path)) => {
            let s = fs::read_to_string(path).map_err(|source| AppError::ReadFile {
                path: path.clone(),
                source,
            })?;
            parse_cmd_input_doc(&s)
        }
        (None, None) => Err(AppError::MissingInputData),
    }
}

type CalculationOutput = CalculationSummary;

#[derive(serde::Deserialize)]
struct CmdInput {
    inputs: Inputs,
    #[serde(default)]
    assumptions: Option<Assumptions>,
}

pub fn print_output(out: &CalculationOutput, args: &Args) -> Result<(), AppError> {
    if args.json {
        let s = serde_json::to_string_pretty(&out)
            .map_err(|source| AppError::SerializeOutput { source })?;
        println!("{}", s);
    } else {
        println!("SP: {:.4}", out.sp);
        println!("SA: {:.4} g/kg", out.sa);
        println!("Density: {:.3} kg/m^3", out.density_kg_per_m3);
        println!("SG 20/20: {:.5}", out.sg_20_20);
        println!("SG 25/25: {:.5}", out.sg_25_25);
    }

    Ok(())
}
