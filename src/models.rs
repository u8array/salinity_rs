use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Assumptions {
    pub temp: f64,
    pub pressure_dbar: f64,
    pub alkalinity: Option<f64>,
    pub alkality: Option<f64>,
    pub assume_borate: bool,
    pub default_f_mg_l: f64,
    pub ref_alk_dkh: Option<f64>,
    pub salinity_norm: f64,
    pub return_components: bool,
    pub borate_fraction: Option<f64>,
    pub alk_mg_per_meq: Option<f64>,
    pub rn_compat: bool,
}

impl Default for Assumptions {
    fn default() -> Self {
        Self {
            temp: 20.0,
            pressure_dbar: 0.0,
            alkalinity: Some(8.0),
            alkality: None,
            assume_borate: true,
            default_f_mg_l: 1.296,
            ref_alk_dkh: Some(8.0),
            salinity_norm: 35.0,
            return_components: false,
            borate_fraction: None,
            alk_mg_per_meq: None,
            rn_compat: false,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Inputs {
    pub na: f64,
    pub ca: f64,
    pub mg: f64,
    pub k: f64,
    pub sr: f64,
    pub br: f64,
    pub cl: Option<f64>,
    pub f: Option<f64>,
    pub s: f64,
    pub b: f64,
    pub t_c: f64,
    pub p_dbar: f64,
    pub alk_dkh: Option<f64>,
    pub assume_borate: bool,
    pub default_f_mg_l: f64,
    pub ref_alk_dkh: Option<f64>,
    pub borate_fraction: Option<f64>,
    pub alk_mg_per_meq: Option<f64>,
    pub return_components: bool,
}

impl Inputs {
    pub fn from_base_with_assumptions(base: &Inputs, ass: &Assumptions) -> Self {
        Self {
            na: base.na,
            ca: base.ca,
            mg: base.mg,
            k: base.k,
            sr: base.sr,
            br: base.br,
            cl: base.cl,
            f: base.f,
            s: base.s,
            b: base.b,
            t_c: ass.temp,
            p_dbar: ass.pressure_dbar,
            alk_dkh: ass.alkalinity,
            assume_borate: ass.assume_borate,
            default_f_mg_l: ass.default_f_mg_l,
            ref_alk_dkh: ass.ref_alk_dkh,
            borate_fraction: ass.borate_fraction,
            alk_mg_per_meq: ass.alk_mg_per_meq,
            return_components: ass.return_components,
        }
    }
}

impl Assumptions {
    pub fn normalized(mut self) -> Self {
        if self.alkalinity.is_none() {
            if let Some(a) = self.alkality {
                self.alkalinity = Some(a);
            }
        }
        if self.rn_compat {
            if self
                .ref_alk_dkh
                .map(|v| (v - crate::chemistry::DEFAULT_REF_ALK_DKH).abs() < f64::EPSILON)
                .unwrap_or(true)
            {
                self.ref_alk_dkh = Some(6.2);
            }
        }
        self
    }
}
