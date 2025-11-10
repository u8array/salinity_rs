// Chemistry constants and helpers
pub const SR_REF: f64 = 35.16504;
pub const DEFAULT_REF_ALK_DKH: f64 = 8.0;

// Molar / mass constants (g/mol)
pub const M_NA: f64 = 22.989_769_28;
pub const M_CA: f64 = 40.078;
pub const M_MG: f64 = 24.305;
pub const M_K: f64 = 39.098_3;
pub const M_SR: f64 = 87.62;
pub const M_BR: f64 = 79.904;
pub const M_CL: f64 = 35.45;
pub const M_F: f64 = 18.998_403_163;
pub const M_S: f64 = 32.065;
pub const M_SO4: f64 = 96.06;
pub const M_B: f64 = 10.81;
pub const M_BORIC: f64 = 61.83; // B(OH)3
pub const M_BORATE: f64 = 60.83; // B(OH)4-

// Reference mmol/kg (standard seawater)
pub const REF_MMOL_CL: f64 = 545.8696;
pub const REF_MMOL_NA: f64 = 468.9674;
pub const REF_MMOL_SO4: f64 = 28.2359;
pub const REF_MMOL_MG: f64 = 52.8116;
pub const REF_MMOL_CA: f64 = 10.2821;
pub const REF_MMOL_K: f64 = 10.2070;
pub const REF_MMOL_BR: f64 = 0.8434;
pub const REF_MMOL_SR: f64 = 0.0906;
pub const REF_MMOL_F: f64 = 0.0680;
pub const REF_MMOL_B: f64 = 0.4160;

// Alkalinity parameters
pub const ALKA_FRAC_HCO3: f64 = 0.89;
pub const ALKA_FRAC_CO3: f64 = 0.10;
pub const ALKA_FRAC_OH: f64 = 0.01;
pub const BORATE_FRACTION_DEFAULT: f64 = 0.20;
pub const DKH_TO_MEQL: f64 = 0.357; // dKH -> meq/L
pub const MG_PER_MEQ_AS_CACO3: f64 = 50.043; // mg/meq as CaCO3

pub const TINY: f64 = 1e-20;
pub const MIN_CL_MG_L: f64 = 0.0;

use crate::models::Inputs;

pub fn sum_ref_gkg() -> f64 {
    let terms = [
        (REF_MMOL_CL, M_CL),
        (REF_MMOL_NA, M_NA),
        (REF_MMOL_SO4, M_SO4),
        (REF_MMOL_MG, M_MG),
        (REF_MMOL_CA, M_CA),
        (REF_MMOL_K, M_K),
        (REF_MMOL_BR, M_BR),
        (REF_MMOL_SR, M_SR),
        (REF_MMOL_F, M_F),
        (REF_MMOL_B, M_B),
    ];
    terms.iter().map(|(mmol, m)| (mmol * m) / 1000.0).sum()
}

pub fn mol_per_l(mg_l: f64, molar_mass_g_mol: f64) -> f64 {
    (mg_l.max(0.0)) / 1000.0 / molar_mass_g_mol.max(TINY)
}

pub fn alk_species_from_dkh(alk_dkh: f64, mg_per_meq: Option<f64>) -> (f64, f64, f64, f64) {
    if alk_dkh <= 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let a_meq_l = alk_dkh * DKH_TO_MEQL;
    let a_eq_l = a_meq_l / 1000.0;

    let a_hco3 = ALKA_FRAC_HCO3 * a_eq_l;
    let a_co3 = ALKA_FRAC_CO3 * a_eq_l;
    let a_oh = ALKA_FRAC_OH * a_eq_l;

    let n_hco3 = a_hco3 / 1.0;
    let n_co3 = a_co3 / 2.0;
    let n_oh = a_oh / 1.0;

    let mg_per_meq_eff = mg_per_meq.unwrap_or(MG_PER_MEQ_AS_CACO3);
    let alk_mass_mg_l = a_meq_l * mg_per_meq_eff;
    (n_hco3, n_co3, n_oh, alk_mass_mg_l)
}

pub fn boron_partition(b_mg_l: f64, borate_fraction: f64) -> (f64, f64) {
    if b_mg_l <= 0.0 {
        return (0.0, 0.0);
    }
    let n_b_total = mol_per_l(b_mg_l, M_B);
    let alpha = borate_fraction.clamp(0.0, 1.0);
    let n_borate = alpha * n_b_total;
    let n_boric = n_b_total - n_borate;
    (n_boric, n_borate)
}

pub fn estimate_cl_mg_l_from_charge_balance(
    inp: &Inputs,
    default_f_mg_l: f64,
    n_borate: f64,
    n_hco3: f64,
    n_co3: f64,
    n_oh: f64,
) -> f64 {
    // Inputs are mg/L; mol_per_l expects mg/L directly.
    let pos = 1.0 * mol_per_l(inp.na, M_NA)
        + 2.0 * mol_per_l(inp.mg, M_MG)
        + 2.0 * mol_per_l(inp.ca, M_CA)
        + 1.0 * mol_per_l(inp.k, M_K)
        + 2.0 * mol_per_l(inp.sr, M_SR);

    // Convert S (mg/L as elemental S) to sulfate mg/L via molar mass ratio.
    let so4_mg_l = (inp.s / M_S) * M_SO4;
    let n_so4 = mol_per_l(so4_mg_l, M_SO4);
    let mut neg = 2.0 * n_so4;
    // Monovalent anions (mg/L)
    neg += mol_per_l(inp.br, M_BR);
    let f_mg_l = inp.f.unwrap_or(default_f_mg_l);
    neg += mol_per_l(f_mg_l, M_F);

    neg += 1.0 * n_borate;
    neg += 1.0 * n_hco3 + 2.0 * n_co3 + 1.0 * n_oh;

    let n_cl = pos - neg;
    let mg_l_cl = (n_cl.max(0.0)) * M_CL * 1000.0;
    mg_l_cl.max(MIN_CL_MG_L)
}

/// Estimate chloride (mg/L) using a blended strategy:
/// 1) Primary: charge balance (as before)
/// 2) Secondary: ion ratio constraints relative to reference seawater composition
///
/// The ratio-based estimate derives n_cl candidates from measured species via
///    n_cl_i = n_i / r_i_ref, where r_i_ref = REF_MMOL_i / REF_MMOL_CL (molar ratios).
/// We then compute a weighted average across available species and blend it with
/// the charge-balance estimate. This improves robustness when Cl is missing.
pub fn estimate_cl_mg_l(
    inp: &Inputs,
    default_f_mg_l: f64,
    n_borate: f64,
    n_hco3: f64,
    n_co3: f64,
    n_oh: f64,
) -> f64 {
    // 1) Charge-balance-based estimate (mol/L)
    let mg_l_charge =
        estimate_cl_mg_l_from_charge_balance(inp, default_f_mg_l, n_borate, n_hco3, n_co3, n_oh);
    let n_cl_charge = (mg_l_charge / 1000.0) / M_CL.max(TINY);

    // 2) Ratio-based candidates (mol/L)
    // Reference molar ratios r_i = REF_MMOL_i / REF_MMOL_CL
    let r_na = REF_MMOL_NA / REF_MMOL_CL;
    let r_mg = REF_MMOL_MG / REF_MMOL_CL;
    let r_ca = REF_MMOL_CA / REF_MMOL_CL;
    let r_k = REF_MMOL_K / REF_MMOL_CL;
    let r_sr = REF_MMOL_SR / REF_MMOL_CL;
    let r_br = REF_MMOL_BR / REF_MMOL_CL;
    let r_so4 = REF_MMOL_SO4 / REF_MMOL_CL;

    // Measured moles per L for species with reliable ratios
    // Inputs are in mg/L – mol_per_l expects mg/L
    let n_na = mol_per_l(inp.na, M_NA);
    let n_mg = mol_per_l(inp.mg, M_MG);
    let n_ca = mol_per_l(inp.ca, M_CA);
    let n_k = mol_per_l(inp.k, M_K);
    let n_sr = mol_per_l(inp.sr, M_SR);
    let n_br = mol_per_l(inp.br, M_BR);
    // Convert S (mg/L as S) to SO4 mg/L by molar mass ratio
    let so4_mg_l = (inp.s / M_S) * M_SO4;
    let n_so4 = mol_per_l(so4_mg_l, M_SO4);

    // Candidate n_cl from each species (ignore invalid/zero)
    let mut sum_w = 0.0;
    let mut sum_w_ncl = 0.0;

    // Weights: proportional to reference molar fraction – Na is strongest, Sr/Br are weaker
    let w_na = REF_MMOL_NA;
    let w_mg = REF_MMOL_MG;
    let w_ca = REF_MMOL_CA;
    let w_k = REF_MMOL_K;
    let w_sr = REF_MMOL_SR;
    let w_br = REF_MMOL_BR;
    let w_so4 = REF_MMOL_SO4;

    let add_candidate = |sum_w: &mut f64, sum_w_ncl: &mut f64, w: f64, n_i: f64, r_i: f64| {
        if w > 0.0 && r_i > 0.0 && n_i > 0.0 {
            let n_cl_i = n_i / r_i;
            *sum_w += w;
            *sum_w_ncl += w * n_cl_i;
        }
    };

    add_candidate(&mut sum_w, &mut sum_w_ncl, w_na, n_na, r_na);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_mg, n_mg, r_mg);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_ca, n_ca, r_ca);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_k, n_k, r_k);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_sr, n_sr, r_sr);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_br, n_br, r_br);
    add_candidate(&mut sum_w, &mut sum_w_ncl, w_so4, n_so4, r_so4);

    let n_cl_ratio = if sum_w > 0.0 {
        (sum_w_ncl / sum_w).max(0.0)
    } else {
        0.0
    };

    // 3) Adaptive blend:
    // If the charge-balance estimate is significantly lower than the ratio-based estimate (underestimation), use the ratio estimate entirely.
    // Otherwise, apply a moderate blend.
    let n_cl_blend = if n_cl_ratio > 0.0 {
        if n_cl_charge < 0.8 * n_cl_ratio {
            n_cl_ratio
        } else {
            let alpha = 0.6;
            alpha * n_cl_charge + (1.0 - alpha) * n_cl_ratio
        }
    } else {
        n_cl_charge
    };

    (n_cl_blend * M_CL * 1000.0).max(MIN_CL_MG_L)
}

pub fn ref_sum_with_boron_species_and_ref_alk(
    ref_alk_dkh: Option<f64>,
    assume_borate: bool,
    borate_fraction_override: Option<f64>,
    alk_mg_per_meq_override: Option<f64>,
) -> f64 {
    let sum_ref_gkg_base = sum_ref_gkg();
    let ref_b_mol_per_kg = REF_MMOL_B / 1000.0;

    let f_borate_ref = if assume_borate {
        borate_fraction_override.unwrap_or(BORATE_FRACTION_DEFAULT)
    } else {
        0.0
    };
    let ref_n_borate = f_borate_ref * ref_b_mol_per_kg;
    let ref_n_boric = ref_b_mol_per_kg - ref_n_borate;

    let ref_b_element_gkg = ref_b_mol_per_kg * M_B;
    let ref_b_species_gkg = ref_n_borate * M_BORATE + ref_n_boric * M_BORIC;

    let ref_alk_gkg = if let Some(val) = ref_alk_dkh {
        if val > 0.0 {
            let (_, _, _, ref_alk_mg_l) = alk_species_from_dkh(val, alk_mg_per_meq_override);
            ref_alk_mg_l / 1000.0
        } else {
            0.0
        }
    } else {
        0.0
    };

    (sum_ref_gkg_base - ref_b_element_gkg) + ref_b_species_gkg + ref_alk_gkg
}

pub fn round_to(x: f64, digits: i32) -> f64 {
    let p = 10f64.powi(digits);
    (x * p).round() / p
}
