use salinity_rs::{
    Assumptions, CalcResult, Inputs, calc_salinity_sp_teos10, rho_from_sp, specific_gravity,
};

fn approx_in_range(v: f64, min: f64, max: f64) {
    assert!((min..=max).contains(&v), "value {v} not in [{min}, {max}]");
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
        t_c: 20.0,
        p_dbar: 0.0,
        alk_dkh: None,
        assume_borate: true,
        default_f_mg_l: 1.296,
        ref_alk_dkh: Some(8.0),
        borate_fraction: None,
        alk_mg_per_meq: None,
        return_components: false,
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
