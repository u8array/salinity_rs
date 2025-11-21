#![allow(clippy::inconsistent_digit_grouping, dead_code)]

#[cfg(feature = "approx_ct")]
use crate::adapters::manual_ct::ct_from_t_manual;
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
/// Falls back to 1025.0 kg/m³ if the `gsw` library returns an error.
pub fn rho(sa: f64, ct: f64, p_dbar: f64) -> f64 {
    gsw_teos10::volume::rho(sa, ct, p_dbar).unwrap_or(1025.0)
}

// --- Reduced TEOS-10 constants (from GSW_TEOS10_CONSTANTS) ---
// These constants mirror those used by the reference TEOS-10 code. We
// duplicate them here to avoid relying on internal symbol names of the
// upstream `gsw` crate; they remain private implementation detail of
// this adapter module.
const GSW_SFAC: f64 = 0.024_882_667_558_461_5;
const GSW_CP0: f64 = 3991.867_957_119_63; // J/(kg·K)
const GSW_T0: f64 = 273.15; // K
const GSW_SSO: f64 = 35.165_04; // Standard Ocean Salinity
const GSW_UPS: f64 = 35.0; // Reference Practical Salinity divisor

fn ct_from_pt(_sa: f64, pt0: f64) -> f64 {
    pt0
}

/// Gibbs second derivative with respect to temperature at p=0.
#[allow(clippy::excessive_precision)]
fn gibbs_pt0_pt0(sa: f64, pt0: f64) -> f64 {
    let x2 = GSW_SFAC * sa;
    let x = x2.sqrt();
    let y = pt0 * 0.025;

    let g03 = -24715.571_866_078
        + y * (4420.447_224_909_672_5
            + y * (-1778.231_237_203_896
                + y * (1160.518_251_685_141_9
                    + y * (-569.531_539_542_516 + y * 128.134_291_524_946_15))));
    let g08 = x2
        * (1760.062_705_994_408
            + x * (-86.132_935_195_608_4
                + x * (-137.114_501_840_898_2
                    + y * (296.200_616_913_752_36
                        + y * (-205.677_092_903_745_63 + 49.939_401_913_901_6 * y)))
                + y * (-60.136_422_517_125 + y * 10.507_207_941_707_34))
            + y * (-1351.605_895_580_406
                + y * (1097.112_537_301_510_9
                    + y * (-433.206_481_750_622_06 + 63.905_091_254_154_904 * y))));
    (g03 + g08) * 0.000_625
}

/// Entropy part at given pressure (TEOS-10 `gsw_entropy_part`).
#[allow(clippy::excessive_precision)]
fn entropy_part(sa: f64, t: f64, p_dbar: f64) -> f64 {
    let x2 = GSW_SFAC * sa;
    let x = x2.sqrt();
    let y = t * 0.025;
    let z = p_dbar * 1e-4;

    let g03 = z
        * (-270.983_805_184_062
            + z * (776.153_611_613_101
                + z * (-196.512_550_881_22
                    + (28.979_652_629_417_5 - 2.132_900_835_183_27 * z) * z)))
        + y * (-24_715.571_866_078
            + z * (2_910.072_908_093_6
                + z * (-1_513.116_771_538_718
                    + z * (546.959_324_647_056
                        + z * (-111.120_812_763_443_6 + 8.688_413_438_343_94 * z)))))
        + y * y
            * (2_210.223_612_454_836_3
                + z * (-2_017.523_349_435_21
                    + z * (1_498.081_172_457_456
                        + z * (-718.635_991_963_235_9
                            + (146.403_755_578_161_6 - 4.989_213_186_267_150_5 * z) * z))))
        + y * y
            * y
            * (-592.743_745_734_632
                + z * (1_591.873_781_627_888
                    + z * (-1_207.261_522_487_504
                        + (608.785_486_935_364 - 105.499_350_893_120_8 * z) * z)))
        + y * y
            * y
            * y
            * (290.129_562_921_285_47
                + z * (-973.091_553_087_975
                    + z * (602.603_274_510_125
                        + z * (-276.361_526_170_076 + 32.409_533_403_861_05 * z))))
        + y * y
            * y
            * y
            * y
            * (-113.906_307_908_503_21
                + y * (21.355_715_254_157_69 - 67.417_568_357_514_34 * z)
                + z * (381.068_361_985_070_96
                    + z * (-133.738_390_284_275_4 + 49.023_632_509_086_724 * z)));

    let g08 = x2
        * (z * (729.116_529_735_046
            + z * (-343.956_902_961_561
                + z * (124.687_671_116_248
                    + z * (-31.656_964_386_073 + 7.046_588_033_154_49 * z))))
            + x * (x
                * (y * (-137.114_501_840_898_2
                    + y * (148.100_308_456_876_18
                        + y * (-68.559_030_967_915_2 + 12.484_850_478_475_4 * y)))
                    - 22.668_355_851_282_9 * z)
                + z * (-175.292_041_186_547
                    + (83.192_392_780_181_9 - 29.483_064_349_429 * z) * z)
                + y * (-86.132_935_195_608_4
                    + z * (766.116_132_004_952
                        + z * (-108.383_452_503_422_4 + 51.279_697_477_982_8 * z))
                    + y * (-30.068_211_258_562_5 - 1_380.959_795_403_770_8 * z
                        + y * (3.502_402_647_235_78 + 938.260_750_445_42 * z))))
            + y * (1_760.062_705_994_408
                + y * (-675.802_947_790_203
                    + y * (365.704_179_100_503_6
                        + y * (-108.301_620_437_655_52 + 12.781_018_250_830_98 * y)
                        + z * (-1_190.914_967_948_748
                            + (298.904_564_555_024 - 145.949_167_600_635_2 * z) * z))
                    + z * (2_082.734_442_399_804_3
                        + z * (-614.668_925_894_709
                            + (340.685_093_521_782 - 33.384_820_297_923_9 * z) * z)))
                + z * (-1_721.528_607_567_954
                    + z * (674.819_060_538_734
                        + z * (-356.629_112_415_276
                            + (88.408_071_661_6 - 15.840_030_944_233_64 * z) * z)))));

    -(g03 + g08) * 0.025
}

fn entropy_part_zerop(sa: f64, pt0: f64) -> f64 {
    entropy_part(sa, pt0, 0.0)
}

#[allow(clippy::excessive_precision)]
fn pt0_from_t(sa: f64, t: f64, p_dbar: f64) -> f64 {
    let s1 = sa / GSW_UPS;
    let mut pt0 = t + p_dbar
        * (8.654_839_133_954_42e-6
            - s1 * 1.416_362_997_448_81e-6
            - p_dbar * 7.382_864_671_357_37e-9
            + t * (-8.382_413_570_396_98e-6
                + s1 * 2.839_333_685_855_34e-8
                + t * 1.778_039_652_186_56e-8
                + p_dbar * 1.711_556_192_082_33e-10));

    let mut dentropy_dt = GSW_CP0 / ((GSW_T0 + pt0) * (1.0 - 0.05 * (1.0 - sa / GSW_SSO)));
    let true_entropy_part = entropy_part(sa, t, p_dbar);

    for _ in 0..2 {
        let pt0_old = pt0;
        let dentropy = entropy_part_zerop(sa, pt0_old) - true_entropy_part;
        pt0 = pt0_old - dentropy / dentropy_dt;
        let pt0m = 0.5 * (pt0 + pt0_old);
        dentropy_dt = -gibbs_pt0_pt0(sa, pt0m);
        pt0 = pt0_old - dentropy / dentropy_dt;
    }
    pt0
}

#[cfg(test)]
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
