use crate::adapters::teos10::{ct_from_t, rho, sa_from_sp};
use crate::chemistry::*;
use crate::models::{Assumptions, Inputs};
use serde::Serialize;

#[derive(Debug)]
pub enum CalcResult {
    Simple(f64),
    Detailed(DetailedResult),
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct Components {
    pub mg_l: Vec<(&'static str, f64)>,
    pub mgkg: Vec<(&'static str, f64)>,
    pub mg_l_sp35: Vec<(&'static str, f64)>,
    pub mgkg_sp35: Vec<(&'static str, f64)>,
    pub norm_factor: f64,
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct DetailedResult {
    pub sp: f64,
    pub rho_kg_m3: f64,
    pub components: Components,
}

#[derive(Serialize, Debug, Clone)]
pub struct CalculationSummary {
    pub sp: f64,
    pub sa: f64,
    pub density_kg_per_m3: f64,
    pub sg_20_20: f64,
    pub sg_25_25: f64,
}

pub fn calc_salinity_sp_iterative(inp: &Inputs, max_iter: usize, tol: f64) -> CalcResult {
    let (n_boric, n_borate) = boron_partition(
        inp.b,
        if inp.assume_borate {
            inp.borate_fraction.unwrap_or(BORATE_FRACTION_DEFAULT)
        } else {
            0.0
        },
    );

    let (n_hco3, n_co3, n_oh, alk_mg_l) =
        alk_species_from_dkh(inp.alk_dkh.unwrap_or(0.0), inp.alk_mg_per_meq);

    let cl_mg_l = inp.cl.filter(|&c| c > 0.0).unwrap_or_else(|| {
        estimate_cl_mg_l_from_charge_balance(inp, n_borate, n_hco3, n_co3, n_oh)
    });

    let f_mg_l = inp.f.unwrap_or(inp.default_f_mg_l);
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
        inp.ref_alk_dkh,
        inp.assume_borate,
        inp.borate_fraction,
        inp.alk_mg_per_meq,
    );

    let mut sp = 35.0;
    let mut sa = sp * (SR_REF / 35.0);
    for _ in 0..max_iter {
        let ct = ct_from_t(sa, inp.t_c, inp.p_dbar);
        let rho_val = rho(sa, ct, inp.p_dbar);
        let kg_per_l = rho_val / 1000.0;
        let sum_meas_gkg: f64 = [
            g_l_na, g_l_ca, g_l_mg, g_l_k, g_l_sr, g_l_br, g_l_f, g_l_so4, g_l_boric, g_l_borate,
            g_l_alk, g_l_cl,
        ]
        .into_iter()
        .sum::<f64>()
            / kg_per_l;

        let sr_new = SR_REF * (sum_meas_gkg / sum_ref_gkg.max(TINY));
        let sp_new = 35.0 * sr_new / SR_REF;
        let sa_new = sr_new;
        if (sp_new - sp).abs() < tol {
            sp = sp_new;
            sa = sa_new;
            break;
        }
        sp = sp_new;
        sa = sa_new;
    }

    if !inp.return_components {
        return CalcResult::Simple(round_to(sp, 4));
    }

    let rho_final = {
        let ct = ct_from_t(sa, inp.t_c, inp.p_dbar);
        rho(sa, ct, inp.p_dbar)
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

pub fn calc_salinity_sp_teos10(
    base_inp: &Inputs,
    ass: &Assumptions,
    max_iter: usize,
    tol: f64,
) -> CalcResult {
    let inp2 = Inputs::from_base_with_assumptions(base_inp, &ass.clone().normalized());
    calc_salinity_sp_iterative(&inp2, max_iter, tol)
}

pub fn rho_from_sp(sp: f64, ass: &Assumptions) -> f64 {
    let sa = sa_from_sp(sp);
    let ct = ct_from_t(sa, ass.temp, ass.pressure_dbar);
    rho(sa, ct, ass.pressure_dbar)
}

pub fn specific_gravity(sp: f64, t_ref: f64, p_ref: f64) -> f64 {
    let sa = sa_from_sp(sp);
    let ct_sw = ct_from_t(sa, t_ref, p_ref);
    let rho_sw = rho(sa, ct_sw, p_ref);
    let ct_pw = ct_from_t(0.0, t_ref, p_ref);
    let rho_pw = rho(0.0, ct_pw, p_ref);
    if rho_pw == 0.0 { 1.0 } else { rho_sw / rho_pw }
}

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
