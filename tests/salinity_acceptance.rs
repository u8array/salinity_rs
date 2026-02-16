use salinity_rs::{
    Assumptions, CalcResult, Inputs, calc_salinity_sp_teos10, rho_from_sp, specific_gravity,
};

fn approx_in_range(v: f64, min: f64, max: f64) {
    assert!((min..=max).contains(&v), "value {v} not in [{min}, {max}]");
}

fn approx_eq(v: f64, expected: f64, tol: f64) {
    assert!(
        (v - expected).abs() <= tol,
        "value {v} differs from expected {expected} by more than {tol}"
    );
}

#[test]
fn calculates_salinity_reasonable_without_explicit_cl() {
    // Same inputs as above but without explicit chloride; the improved estimator
    // should still land in a reasonable seawater salinity range.
    let inputs = Inputs {
        na: 11_980.0,
        ca: 357.0,
        mg: 1_246.0,
        k: 464.0,
        sr: 6.96,
        br: 73.2,
        cl: None,
        f: Some(1.14),
        s: 814.0,
        b: 5.57,
        alk_dkh: None,
    };

    let ass = Assumptions {
        ..Default::default()
    };

    let sp = match calc_salinity_sp_teos10(&inputs, &ass, 30, 1e-8) {
        CalcResult::Simple(v) => v,
        CalcResult::Detailed(d) => d.sp,
    };

    // With the blended estimator, SP should be within typical seawater band
    approx_in_range(sp, 30.0, 40.0);
}

#[test]
fn calculates_salinity_for_sample_inputs_within_reasonable_bounds() {
    let inputs = Inputs {
        na: 11_980.0,
        ca: 357.0,
        mg: 1_246.0,
        k: 464.0,
        sr: 6.96,
        br: 73.2,
        cl: Some(19_570.0),
        f: Some(1.14),
        s: 814.0,
        b: 5.57,
        alk_dkh: None,
    };

    let ass = Assumptions {
        ..Default::default()
    };

    let sp = match calc_salinity_sp_teos10(&inputs, &ass, 30, 1e-8) {
        CalcResult::Simple(v) => v,
        CalcResult::Detailed(d) => d.sp,
    };

    // Plausibility checks instead of hard snapshots (more robust against small numerical drifts)
    approx_in_range(sp, 30.0, 40.0); // typical seawater salinity

    let rho = rho_from_sp(sp, &ass);
    approx_in_range(rho, 1000.0, 1060.0); // density (kg/m^3)

    let sg_20 = specific_gravity(sp, 20.0, 0.0);
    approx_in_range(sg_20, 0.98, 1.10);
}

#[test]
fn salinity_norm_changes_component_normalization_factor() {
    let inputs = Inputs {
        na: 11_980.0,
        ca: 357.0,
        mg: 1_246.0,
        k: 464.0,
        sr: 6.96,
        br: 73.2,
        cl: Some(19_570.0),
        f: Some(1.14),
        s: 814.0,
        b: 5.57,
        alk_dkh: None,
    };

    let ass35 = Assumptions {
        return_components: true,
        salinity_norm: 35.0,
        ..Default::default()
    };
    let ass50 = Assumptions {
        return_components: true,
        salinity_norm: 50.0,
        ..Default::default()
    };

    let factor35 = match calc_salinity_sp_teos10(&inputs, &ass35, 30, 1e-8) {
        CalcResult::Detailed(d) => d.components.norm_factor,
        CalcResult::Simple(_) => panic!("expected detailed output"),
    };
    let factor50 = match calc_salinity_sp_teos10(&inputs, &ass50, 30, 1e-8) {
        CalcResult::Detailed(d) => d.components.norm_factor,
        CalcResult::Simple(_) => panic!("expected detailed output"),
    };

    assert!((factor50 / factor35 - (50.0 / 35.0)).abs() < 1e-12);
}

#[test]
fn summary_matches_reference_values_for_known_sample() {
    let inputs = Inputs {
        na: 11_980.0,
        ca: 357.0,
        mg: 1_246.0,
        k: 464.0,
        sr: 6.96,
        br: 73.2,
        cl: Some(19_570.0),
        f: Some(1.14),
        s: 814.0,
        b: 5.57,
        alk_dkh: None,
    };
    let ass = Assumptions {
        ..Default::default()
    };

    let summary = salinity_rs::compute_summary(&inputs, &ass);

    approx_eq(summary.sp, 35.2417, 1e-4);
    approx_eq(summary.sa, 35.407_879_719_085_71, 1e-9);
    approx_eq(summary.density_kg_per_m3, 1_024.945_755_987_009_5, 1e-9);
    approx_eq(summary.sg_20_20, 1.026_579_948_021_515_4, 1e-12);
    approx_eq(summary.sg_25_25, 1.026_237_065_651_817_4, 1e-12);
}
