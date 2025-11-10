use salinity_rs::{Assumptions, Inputs, chemistry::*};

#[test]
fn estimates_chloride_close_to_reference_when_missing() {
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
        alk_dkh: Some(8.0),
    };
    let ass = Assumptions {
        ..Default::default()
    };

    // Partition boron and alkalinity species as solver would
    let (_n_boric, n_borate) = boron_partition(
        inputs.b,
        ass.borate_fraction.unwrap_or(BORATE_FRACTION_DEFAULT),
    );
    let (n_hco3, n_co3, n_oh, _alk_mg_l) =
        alk_species_from_dkh(ass.alkalinity.unwrap_or(8.0), ass.alk_mg_per_meq);

    let cl_mg_l = estimate_cl_mg_l(&inputs, ass.default_f_mg_l, n_borate, n_hco3, n_co3, n_oh);

    // Seawater chloride (mg/L) at SP~35 typically ~19,000–21,000 mg/L depending on composition/density.
    assert!(
        cl_mg_l > 18000.0 && cl_mg_l < 23000.0,
        "estimated Cl mg/L unexpected: {cl_mg_l}"
    );
}
