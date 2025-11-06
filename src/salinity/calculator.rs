use crate::adapters::teos10::{ct_from_t, rho, sa_from_sp};
use crate::chemistry::*;
use crate::models::{Assumptions, Inputs};
use serde::Serialize;

/// Result of a salinity calculation.
///
/// This enum expresses the two possible return shapes from the calculation
/// functions in this module:
///
/// - `Simple(f64)` — a single practical salinity (SP) value.
/// - `Detailed(DetailedResult)` — a richer result with density and per-ion
///   component breakdowns.
#[derive(Debug)]
pub enum CalcResult {
    /// A single SP value (practical salinity).
    Simple(f64),
    /// A detailed result including density and component concentrations.
    Detailed(DetailedResult),
}

/// Concentration tables for ionic and related components.
///
/// All concentration vectors hold tuples of `(name, value)` where the value
/// units are indicated by the field name:
/// - `mg_l`: milligrams per litre (mg/L) at the sample density
/// - `mgkg`: milligrams per kilogram (mg/kg) at the sample density
/// - `mg_l_sp35`: mg/L normalized to SP = 35
/// - `mgkg_sp35`: mg/kg normalized to SP = 35
///
/// `norm_factor` is the multiplicative factor used to normalize component
/// values to SP = 35.
#[derive(Debug)]
pub struct Components {
    pub mg_l: Vec<(&'static str, f64)>,
    pub mgkg: Vec<(&'static str, f64)>,
    pub mg_l_sp35: Vec<(&'static str, f64)>,
    pub mgkg_sp35: Vec<(&'static str, f64)>,
    pub norm_factor: f64,
}

/// A detailed result returned when the caller requests component output.
///
/// - `sp`: practical salinity (rounded in the caller before return)
/// - `rho_kg_m3`: in-situ density in kg/m^3
/// - `components`: per-ion concentration tables and normalization factor
#[derive(Debug)]
pub struct DetailedResult {
    pub sp: f64,
    pub rho_kg_m3: f64,
    pub components: Components,
}

/// Lightweight summary returned for higher-level callers (e.g. UI or API).
///
/// Fields:
/// - `sp`: practical salinity
/// - `sa`: absolute salinity (g/kg)
/// - `density_kg_per_m3`: in-situ density at the sample conditions
/// - `sg_20_20`: specific gravity at 20°C/20°C reference (unitless)
/// - `sg_25_25`: specific gravity at 25°C/25°C reference (unitless)
#[derive(Serialize, Debug, Clone)]
pub struct CalculationSummary {
    pub sp: f64,
    pub sa: f64,
    pub density_kg_per_m3: f64,
    pub sg_20_20: f64,
    pub sg_25_25: f64,
}

/// Compute practical salinity (SP) iteratively from an `Inputs` structure.
///
/// The function builds a mass budget from the provided ion concentrations and
/// iteratively adjusts SP (and therefore the salinity ratio) so that the
/// measured sum of major ions matches the reference sum for the chosen
/// reference composition. When `ass.return_components` is `false` a compact
/// `CalcResult::Simple(sp)` is returned; otherwise a `CalcResult::Detailed`
/// with density and per-ion breakdowns is returned.
///
/// Parameters:
/// - `inp`: computation inputs (concentrations, temperature, pressure, flags).
/// - `max_iter`: maximum number of iterations of the SP update loop.
/// - `tol`: convergence tolerance applied to SP changes.
///
/// Returns: `CalcResult` (either `Simple(f64)` or `Detailed(DetailedResult)`).
pub fn calc_salinity_sp_iterative(
    inp: &Inputs,
    ass: &Assumptions,
    max_iter: usize,
    tol: f64,
) -> CalcResult {
    // Partition boron between boric acid and borate based on assumptions.
    let (n_boric, n_borate) = boron_partition(
        inp.b,
        if ass.assume_borate {
            ass.borate_fraction.unwrap_or(BORATE_FRACTION_DEFAULT)
        } else {
            0.0
        },
    );

    // Convert alkalinity (DKH or mg per meq) into species and total alkalinity
    // in mg/L for the mass-balance. The tuple contains derived species and
    // the equivalent alkalinity in mg/L used directly below.
    let alk_dkh_eff = inp.alk_dkh.or(ass.alkalinity).unwrap_or(0.0);
    let (n_hco3, n_co3, n_oh, alk_mg_l) = alk_species_from_dkh(alk_dkh_eff, ass.alk_mg_per_meq);

    // Chloride: use provided value if positive, otherwise estimate by charge
    // balance using the derived anion/cation species above.
    let cl_mg_l = inp.cl.filter(|&c| c > 0.0).unwrap_or_else(|| {
        estimate_cl_mg_l_from_charge_balance(inp, ass.default_f_mg_l, n_borate, n_hco3, n_co3, n_oh)
    });

    let f_mg_l = inp.f.unwrap_or(ass.default_f_mg_l);
    let g_l_na = inp.na.max(0.0) / 1000.0;
    let g_l_ca = inp.ca.max(0.0) / 1000.0;
    let g_l_mg = inp.mg.max(0.0) / 1000.0;
    let g_l_k = inp.k.max(0.0) / 1000.0;
    let g_l_sr = inp.sr.max(0.0) / 1000.0;
    let g_l_br = inp.br.max(0.0) / 1000.0;
    let g_l_f = f_mg_l.max(0.0) / 1000.0;
    let g_l_so4 = ((inp.s / M_S) * M_SO4 / 1000.0).max(0.0);
    let g_l_boric = n_boric * M_BORIC;
    let g_l_borate = n_borate * M_BORATE;
    let g_l_alk = alk_mg_l / 1000.0;
    let g_l_cl = cl_mg_l.max(0.0) / 1000.0;

    let sum_ref_gkg = ref_sum_with_boron_species_and_ref_alk(
        ass.ref_alk_dkh,
        ass.assume_borate,
        ass.borate_fraction,
        ass.alk_mg_per_meq,
    );

    // Start the iteration from a nominal SP = 35 and update until convergence.
    // The iteration adjusts the salinity ratio so that the measured sum of
    // dissolved species (converted to g/kg) matches the reference composition
    // scaled by that salinity ratio.
    let mut sp = 35.0;
    let mut sa = sp * (SR_REF / 35.0);
    for _ in 0..max_iter {
        // Compute conservative temperature and density at current SA.
        let ct = ct_from_t(sa, ass.temp, ass.pressure_dbar);
        let rho_val = rho(sa, ct, ass.pressure_dbar);
        let kg_per_l = rho_val / 1000.0;

        // Sum the provided mass contributions (g/L) and convert to g/kg by
        // dividing by the in-situ kg/L.
        let sum_meas_gkg: f64 = [
            g_l_na, g_l_ca, g_l_mg, g_l_k, g_l_sr, g_l_br, g_l_f, g_l_so4, g_l_boric, g_l_borate,
            g_l_alk, g_l_cl,
        ]
        .into_iter()
        .sum::<f64>()
            / kg_per_l;

        // New salinity ratio (sr) is set by matching the measured sum to the
        // reference sum (with protection against division by tiny values).
        let sr_new = SR_REF * (sum_meas_gkg / sum_ref_gkg.max(TINY));
        let sp_new = 35.0 * sr_new / SR_REF;
        let sa_new = sr_new;

        // Check for convergence on the practical salinity (SP).
        if (sp_new - sp).abs() < tol {
            sp = sp_new;
            sa = sa_new;
            break;
        }
        sp = sp_new;
        sa = sa_new;
    }

    // If the caller did not request component output, return a compact value.
    if !ass.return_components {
        return CalcResult::Simple(round_to(sp, 4));
    }

    // Recompute final density at the converged SA for output.
    let rho_final = {
        let ct = ct_from_t(sa, ass.temp, ass.pressure_dbar);
        rho(sa, ct, ass.pressure_dbar)
    };
    let kg_per_l = rho_final / 1000.0;

    let mg_l_table = vec![
        ("Na+", g_l_na * 1000.0),
        ("Ca2+", g_l_ca * 1000.0),
        ("Mg2+", g_l_mg * 1000.0),
        ("K+", g_l_k * 1000.0),
        ("Sr2+", g_l_sr * 1000.0),
        ("Br-", g_l_br * 1000.0),
        ("SO4^2-", g_l_so4 * 1000.0),
        ("F-", g_l_f * 1000.0),
        ("Alk.", g_l_alk * 1000.0),
        ("B(OH)3", g_l_boric * 1000.0),
        ("B(OH)4-", g_l_borate * 1000.0),
        ("Cl-", g_l_cl * 1000.0),
    ];

    let mgkg_table: Vec<(&str, f64)> = mg_l_table
        .iter()
        .map(|(k, v)| (*k, *v / kg_per_l))
        .collect();

    let norm_factor = 35.0 / sp.max(TINY);
    let mg_l_norm: Vec<(&str, f64)> = mg_l_table
        .iter()
        .map(|(k, v)| (*k, *v * norm_factor))
        .collect();
    let mgkg_norm: Vec<(&str, f64)> = mgkg_table
        .iter()
        .map(|(k, v)| (*k, *v * norm_factor))
        .collect();

    CalcResult::Detailed(DetailedResult {
        sp: round_to(sp, 4),
        rho_kg_m3: rho_final,
        components: Components {
            mg_l: mg_l_table,
            mgkg: mgkg_table,
            mg_l_sp35: mg_l_norm,
            mgkg_sp35: mgkg_norm,
            norm_factor,
        },
    })
}

/// Compute practical salinity (SP) using TEOS-10 assumptions.
///
/// This is a thin wrapper that combines the provided `base_inp` with the
/// given `Assumptions` (temperature/pressure and other normalized flags)
/// and calls `calc_salinity_sp_iterative` to perform the actual mass-balance
/// based iterative calculation.
///
/// Parameters:
/// - `base_inp`: input concentrations and options (see `Inputs`).
/// - `ass`: calculation `Assumptions` such as temperature and pressure.
/// - `max_iter`: maximum number of iterations for the underlying solver.
/// - `tol`: convergence tolerance applied to SP changes.
///
/// Returns a `CalcResult`, either `Simple(sp)` or `Detailed(...)` depending
/// on the `Inputs` flags.
pub fn calc_salinity_sp_teos10(
    base_inp: &Inputs,
    ass: &Assumptions,
    max_iter: usize,
    tol: f64,
) -> CalcResult {
    let ass_norm = ass.clone().normalized();
    calc_salinity_sp_iterative(base_inp, &ass_norm, max_iter, tol)
}

/// Compute in-situ density (kg/m³) from Practical Salinity (SP).
///
/// The function converts `sp` to absolute salinity (SA) via `sa_from_sp`,
/// computes conservative temperature from the provided `Assumptions` and then
/// calls the TEOS-10 `rho` routine. Units:
/// - `sp` is unitless (practical salinity)
/// - `ass.temp` is °C, `ass.pressure_dbar` is in dbar
/// - return value is density in kg/m³
///
/// This helper is convenient when callers only have SP and a set of
/// environmental assumptions.
pub fn rho_from_sp(sp: f64, ass: &Assumptions) -> f64 {
    let sa = sa_from_sp(sp);
    let ct = ct_from_t(sa, ass.temp, ass.pressure_dbar);
    rho(sa, ct, ass.pressure_dbar)
}

/// Compute the specific gravity of a seawater sample relative to pure water.
///
/// The specific gravity is defined here as the ratio rho(sw) / rho(pw) where
/// `rho(sw)` is the density of seawater at the specified SP, temperature and
/// pressure, and `rho(pw)` is the density of pure water at the same
/// temperature and pressure. If the pure-water density routine returns zero
/// (unexpected), the function conservatively returns `1.0`.
///
/// Parameters:
/// - `sp`: practical salinity (unitless)
/// - `t_ref`: reference temperature in °C
/// - `p_ref`: reference pressure in dbar
///
/// Returns the dimensionless specific gravity (unitless).
pub fn specific_gravity(sp: f64, t_ref: f64, p_ref: f64) -> f64 {
    let sa = sa_from_sp(sp);
    let ct_sw = ct_from_t(sa, t_ref, p_ref);
    let rho_sw = rho(sa, ct_sw, p_ref);
    let ct_pw = ct_from_t(0.0, t_ref, p_ref);
    let rho_pw = rho(0.0, ct_pw, p_ref);
    if rho_pw == 0.0 { 1.0 } else { rho_sw / rho_pw }
}

/// Compute a compact `CalculationSummary` for the given inputs.
///
/// This convenience function runs the TEOS-10 based SP solver and returns a
/// small summary useful for UI or API responses. The returned `CalculationSummary`
/// contains both salinity (SP, SA), the in-situ density (kg/m³) and two
/// reference specific gravities at 20°C and 25°C (both at 0 dbar).
///
/// Notes:
/// - The function uses `calc_salinity_sp_teos10` with conservative defaults
///   for iteration (30 iterations and 1e-8 tolerance).
/// - All units follow the crate convention: density in kg/m³, SA in g/kg,
///   and specific gravities are unitless ratios.
pub fn compute_summary(inputs: &Inputs, assumptions: &Assumptions) -> CalculationSummary {
    let sp = match calc_salinity_sp_teos10(inputs, assumptions, 30, 1e-8) {
        CalcResult::Simple(v) => v,
        CalcResult::Detailed(d) => d.sp,
    };
    let sa = sa_from_sp(sp);
    let rho_val = rho_from_sp(sp, assumptions);
    let sg_20 = specific_gravity(sp, 20.0, 0.0);
    let sg_25 = specific_gravity(sp, 25.0, 0.0);

    CalculationSummary {
        sp,
        sa,
        density_kg_per_m3: rho_val,
        sg_20_20: sg_20,
        sg_25_25: sg_25,
    }
}
