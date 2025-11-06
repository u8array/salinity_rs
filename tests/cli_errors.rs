use predicates::prelude::*;

#[test]
fn cli_fails_without_any_input() {
    let mut cmd = assert_cmd::cargo::cargo_bin_cmd!("salinity_rs");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Missing input data"));
}

#[test]
fn cli_works_without_assumptions_with_inputs_json() {
    let mut cmd = assert_cmd::cargo::cargo_bin_cmd!("salinity_rs");
    let inputs = serde_json::json!({
        "na": 11980.0,
        "ca": 357.0,
        "mg": 1246.0,
        "k": 464.0,
        "sr": 6.96,
        "br": 73.2,
        "cl": 19570.0,
        "f": 1.14,
        "s": 814.0,
        "b": 5.57,
        "t_c": 20.0,
        "p_dbar": 0.0,
        "alk_dkh": null,
        "assume_borate": true,
        "default_f_mg_l": 1.296,
        "ref_alk_dkh": 8.0,
        "borate_fraction": null,
        "alk_mg_per_meq": null,
        "return_components": false,
    })
    .to_string();

    cmd.arg("--json").arg("--inputs-json").arg(inputs);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("\"sp\""));
}

#[test]
fn cli_works_without_assumptions_in_stdin_input_document() {
    let mut cmd = assert_cmd::cargo::cargo_bin_cmd!("salinity_rs");

    let doc = serde_json::json!({
        "inputs": {
            "na": 11980.0,
            "ca": 357.0,
            "mg": 1246.0,
            "k": 464.0,
            "sr": 6.96,
            "br": 73.2,
            "cl": 19570.0,
            "f": 1.14,
            "s": 814.0,
            "b": 5.57,
            "t_c": 20.0,
            "p_dbar": 0.0,
            "alk_dkh": null,
            "assume_borate": true,
            "default_f_mg_l": 1.296,
            "ref_alk_dkh": 8.0,
            "borate_fraction": null,
            "alk_mg_per_meq": null,
            "return_components": false
        }
    })
    .to_string();

    cmd.arg("--json").arg("--input").arg("-").write_stdin(doc);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("\"sp\""));
}

#[test]
fn cli_reports_invalid_json_for_inputs_json() {
    let mut cmd = assert_cmd::cargo::cargo_bin_cmd!("salinity_rs");
    cmd.arg("--inputs-json").arg("{not valid json}");

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Invalid JSON for --inputs-json"));
}

#[test]
fn cli_reports_invalid_json_in_file() {
    use std::fs::File;
    use std::io::Write as _;
    use tempfile::tempdir;

    let dir = tempdir().unwrap();
    let file_path = dir.path().join("bad.json");
    let mut f = File::create(&file_path).unwrap();
    writeln!(f, "this is not json").unwrap();

    let mut cmd = assert_cmd::cargo::cargo_bin_cmd!("salinity_rs");
    cmd.arg("--input").arg(file_path);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Invalid JSON in input document"));
}
