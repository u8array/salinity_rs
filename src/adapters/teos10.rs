#[cfg(feature = "approx_ct")]
use crate::adapters::manual_ct::ct_from_t_manual;
#[cfg(not(feature = "approx_ct"))]
use crate::adapters::teos10_reduced::{ct_from_pt, pt0_from_t};
use gsw as gsw_teos10;

/// Absolute/Reference Salinity from Practical Salinity.
/// Note: This returns TEOS-10 Reference Salinity (SR) from SP and is used
/// as an approximation for Absolute Salinity (SA). For standard seawater
/// composition SR ≈ SA; use location-based SA conversions if available.
pub fn sa_from_sp(sp: f64) -> f64 {
    gsw_teos10::conversions::sr_from_sp(sp)
}

/// Computes Conservative Temperature (CT) from in-situ temperature `t` and Absolute Salinity `sa`.
///
/// # Arguments
/// * `sa` - Absolute Salinity [g/kg]
/// * `temp` - In-situ temperature [°C]
/// * `p_dbar` - Pressure [dbar]
///
/// # Returns
/// * Conservative Temperature [°C]
///
/// # Features
/// * If the feature `reduced_ct` is enabled, uses a fast manual approximation for CT.
/// * Otherwise, calculates a reduced TEOS-10 potential temperature (`pt0`) and applies a simplified identity conversion (CT ≈ PT0).
pub fn ct_from_t(sa: f64, temp: f64, p_dbar: f64) -> f64 {
    #[cfg(feature = "approx_ct")]
    {
        ct_from_t_manual(sa, temp, p_dbar)
    }

    #[cfg(not(feature = "approx_ct"))]
    {
        let pt0 = pt0_from_t(sa, temp, p_dbar);
        ct_from_pt(sa, pt0)
    }
}

/// In-situ density ρ from SA, CT and p (TEOS-10, 75-term polynomial).
/// Returns `NaN` if the `gsw` library reports an error.
pub fn rho(sa: f64, ct: f64, p_dbar: f64) -> f64 {
    gsw_teos10::volume::rho(sa, ct, p_dbar).unwrap_or(f64::NAN)
}

#[cfg(all(test, not(feature = "approx_ct")))]
mod tests {
    use super::*;

    #[test]
    fn pt0_iteration_converges_reasonably() {
        let sa = 35.2;
        let t = 12.0;
        let p = 500.0;
        let pt0 = pt0_from_t(sa, t, p);
        assert!(pt0 < t);
        assert!((t - pt0) < 3.0);
    }

    #[test]
    fn ct_identity_for_now() {
        let sa = 35.0;
        let t = 20.0;
        let p = 10.0;
        let ct = ct_from_t(sa, t, p);
        assert!(ct <= t);
    }
}
