
# salinity_rs

![Build Status](https://github.com/u8array/salinity_rs/actions/workflows/rust.yml/badge.svg)
[![no_std supported](https://img.shields.io/badge/no__std-supported-blue.svg)](Cargo.toml)
<!-- [![crates.io](https://img.shields.io/crates/v/salinity_rs.svg)](https://crates.io/crates/salinity_rs)
[![docs.rs](https://docs.rs/salinity_rs/badge.svg)](https://docs.rs/salinity_rs) -->
[![license](https://img.shields.io/github/license/u8array/salinity_rs.svg)](LICENSE)

Rust library (with optional CLI) to estimate **Practical Salinity** `SP`, **Absolute Salinity** `SA`, **density** `ρ`, and specific gravities from macro chemical analyses of seawater or reef aquaria.

Calculations follow TEOS‑10 conventions and couple charge/mass balances with density via a simple fixed‑point iteration. The library exposes both high‑level helpers and low‑level building blocks; the CLI prints a concise summary (SP, SA, ρ, SG).

## Features

- Compute `SP` and `SA` from elemental/ionic analyses
- TEOS-10 density `ρ(SA,CT,p)` and specific gravity at a reference temperature
- Self-consistent conversion between mg/L and mg/kg via `ρ`
- Chloride estimation from electroneutrality with adaptive reference‑ion ratio blending when Cl⁻ is missing
- Configurable assumptions: temperature `T`, pressure `p`, alkalinity, borate fraction
- Library support for component tables (mg/L, mg/kg, SP=35 normalization). Note: current CLI prints summary only.

## Install and use as a library

Add the crate to your project. Until published on crates.io, you can depend on the Git repository. Disable default features to avoid pulling the optional CLI dependency:

```toml
[dependencies]
salinity_rs = { git = "https://github.com/u8array/salinity_rs", default-features = false }
```

Quick start (Rust):

```rust
use salinity_rs::{compute_summary, Inputs, Assumptions};

fn main() {
  // concentrations in mg/L; temperature in °C; pressure in dbar
  let inputs = Inputs {
    na: 10780.0, ca: 420.0, mg: 1290.0, k: 400.0, sr: 8.0,
    br: 65.0, s: 900.0, b: 4.4,
    cl: None,    // let the model estimate Cl⁻ from electroneutrality
    f: None,     // fall back to default F⁻ if not provided
    alk_dkh: Some(8.0),
  };

  // Environmental and reference assumptions
  let ass = Assumptions { temp: 20.0, pressure_dbar: 0.0, ..Default::default() };
  let out = compute_summary(&inputs, &ass);
  println!(
    "SP={:.4} SA={:.4} g/kg  ρ={:.3} kg/m³  SG20/20={:.5} SG25/25={:.5}",
    out.sp, out.sa, out.density_kg_per_m3, out.sg_20_20, out.sg_25_25
  );
}
```

Public API highlights:

- Convenience API: `compute_summary(inputs, assumptions)` → SP, SA, ρ, SG(20/20), SG(25/25)
- Solver API: `calc_salinity_sp_teos10`, `calc_salinity_sp_iterative(&Inputs, &Assumptions, max_iter, tol)`, `rho_from_sp`, `specific_gravity`, `sa_from_sp`
- Types: `Inputs`, `Assumptions`, `CalcResult`, `DetailedResult`, `Components`

Minimum supported Rust: a recent stable with Edition 2024 support.

## CLI usage (optional feature)

The CLI is behind the `cli` feature. Build and run locally:

```bash
cargo run --features cli -- --help
```

The binary accepts JSON inputs for measurements and assumptions:

- `--inputs-json <JSON>`: Inline JSON for input values (shape of `Inputs`).
- `--assumptions-json <JSON>`: Optional, adds/overrides assumptions (shape of `Assumptions`).
- `--input <FILE>`: Read a file containing an object with `inputs` and optional `assumptions`. Use `-` for stdin.
- `--json`: Output machine‑readable JSON.

JSON fields (excerpt):

- Inputs (mg/L unless noted): `na, ca, mg, k, sr, br, s, b, cl` (optional; omit or set `null` for auto‑estimate), `f` (optional), `alk_dkh` (dKH, optional)
- Assumptions (conditions and options): `temp` (°C, default 20), `pressure_dbar` (dbar, default 0), `alkalinity` (dKH, optional; used if `inputs.alk_dkh` is missing), `assume_borate` (default true), `borate_fraction`, `ref_alk_dkh` (default 8.0), `alk_mg_per_meq`, `default_f_mg_l` (default 1.296), `return_components` (default false)

## Output example

```text
SP: 34.9800
SA: 35.1600 g/kg
Density: 1024.600 kg/m^3
SG 20/20: 1.02600
SG 25/25: 1.02480
```

## Assumptions and limits

- $\delta SA_{\text{composition}}$ neglected; acceptable for near-NSW compositions.
- Fixed alkalinity fractions are pH-independent; accuracy drops for unusual pH/CO₂.
- Cl⁻ by electroneutrality is sensitive to input errors.
- Optional default F⁻ used if missing.
- Iteration enforces consistency between volume-based inputs and mass-based reference.

## Feature flags

- `cli` — enables the command‑line interface and pulls in the optional `clap` dependency. Not needed for library use.
- `approx_ct` — enables a shallow‑water approximation for CT from t (CT≈θ via small adiabatic lapse). For aquaria/near‑surface use only; otherwise leave disabled (default).
- `std` — enables use of the Rust standard library; when disabled the crate can be built with `no_std`.

## Quick start

Example with concentrations in mg/L, `t = 20 °C`, `p = 0 dbar`, and total alkalinity `8 dKH`.

Two invocation methods are supported:

1. Inline JSON (on the command line)

```bash
salinity_rs --inputs-json '{
  "na":10780, "mg":1290, "ca":420, "k":400, "sr":8,
  "br":65, "s":900, "b":4.4,
  "alk_dkh":8
}' --assumptions-json '{
  "temp":20, "pressure_dbar":0, "return_components":true
}'
```

1. File-based (file contains an object with `inputs` and optional `assumptions`)

`inputs.json`

```json
{
  "inputs": {
    "na":10780, "mg":1290, "ca":420, "k":400, "sr":8,
    "br":65, "s":900, "b":4.4,
    "alk_dkh":8
  },
  "assumptions": {
    "temp": 20, "pressure_dbar": 0, "return_components": true
  }
}
```

Run:

```bash
salinity_rs --input inputs.json
```

## Background

### Practical vs Absolute Salinity

- **SP**: dimensionless, historically defined from conductivity relative to standard seawater at 15 °C and 1 bar.
- **SA**: g/kg of dissolved material. Primary TEOS-10 salinity variable.  
  For near-standard composition:

```math
SA \approx SR = SP\cdot\frac{SR_{\mathrm{REF}}}{35},\qquad SR_{\mathrm{REF}}=35.16504\ \mathrm{g\ kg^{-1}}.
```

  Composition anomalies would add $\delta SA_{\text{composition}}$. This tool assumes $\delta SA_{\text{composition}}\approx 0$.

### Units and basic conversions

For species $i$ with input concentration $c_i$ in mg/L and molar mass $M_i$ in g/mol:

```math
n_i\ [\mathrm{mol\ L^{-1}}] = \frac{\max(c_i,0)}{1000}\cdot\frac{1}{M_i}.
```

### Alkalinity split (approximate)

Total alkalinity in dKH:

```math
A_{\mathrm{meq/L}} = \mathrm{dKH}\cdot 0.357.
```

Split into species by fixed fractions:

```math
(A_{\mathrm{HCO_3^-}}, A_{\mathrm{CO_3^{2-}}}, A_{\mathrm{OH^-}}) = (0.89,\,0.10,\,0.01)\cdot A_{\mathrm{meq/L}},
```

with stoichiometry

```math
n_{\mathrm{HCO_3^-}}=A_{\mathrm{HCO_3^-}},\quad
n_{\mathrm{CO_3^{2-}}}=\frac{A_{\mathrm{CO_3^{2-}}}}{2},\quad
n_{\mathrm{OH^-}}=A_{\mathrm{OH^-}}.
```

Mass equivalent for reporting uses a configurable $\mathrm{mg/meq}$ (default 50.043 mg/meq as CaCO₃).  
Note: exact speciation is pH and DIC dependent; fixed fractions are a robust aquarium approximation.

### Boron partition

Total B split by borate fraction $\alpha_B\in[0,1]$:

```math
n_{\mathrm{B(OH)_4^-}}=\alpha_B\,n_B,\qquad
n_{\mathrm{B(OH)_3}}=(1-\alpha_B)\,n_B.
```

This affects charge and mass because species have different molar masses.

### Electroneutrality and Cl⁻ estimation

If Cl⁻ is not provided, solve from

```math
\sum_i z_i n_i = 0.
```

Cations: $\mathrm{Na^+}, \mathrm{Mg^{2+}}, \mathrm{Ca^{2+}}, \mathrm{K^+}, \mathrm{Sr^{2+}}$.  
Anions: $\mathrm{SO_4^{2-}}, \mathrm{Br^-}, \mathrm{F^-}, \mathrm{B(OH)_4^-}, \mathrm{HCO_3^-}, \mathrm{CO_3^{2-}}, \mathrm{OH^-}$.  
Assign the residual negative charge to $\mathrm{Cl^-}$ and clamp at zero if needed.

#### Adaptive chloride estimation (charge balance + ratio blending)

When $\mathrm{Cl^-}$ is not measured, the library complements the charge-balance estimate with a ratio-based estimate derived from the reference composition. Let the reference molar ratios be

```math
r_i = \frac{\mathrm{REF\_MMOL}_i}{\mathrm{REF\_MMOL}_{\mathrm{Cl}}}.
```

From measured mass concentrations $c_i$ (mg/L) and molar masses $M_i$ (g/mol), obtain $n_i$ (mol/L) via

```math
n_i = \frac{\max(c_i,0)}{1000\,M_i}.
```

For each species with a sufficiently stable ratio in near‑NSW (Na⁺, Mg²⁺, Ca²⁺, K⁺, Sr²⁺, Br⁻, and SO₄²⁻), form a candidate

```math
n_{\mathrm{Cl},i} = \frac{n_i}{r_i}.
```

Use reference mmol/kg as weights $w_i$ (Na strongest influence) to compute a weighted mean

```math
n_{\mathrm{Cl,ratio}} = \frac{\sum_i w_i\,n_{\mathrm{Cl},i}}{\sum_i w_i},\qquad w_i=\mathrm{REF\_MMOL}_i.
```

Derive $n_{\mathrm{SO_4}}$ from sulfur reported as elemental S using the molar mass ratio, i.e., $c_{\mathrm{SO_4}} = c_S\,M_{\mathrm{SO_4}}/M_S$ prior to conversion to moles. Fluoride is handled in charge balance (with a configurable default if missing) but is not used in the ratio estimate.

Finally, blend the charge‑balance and ratio estimates adaptively:

```math
n_{\mathrm{Cl,blend}} =
\begin{cases}
 n_{\mathrm{Cl,ratio}}, & n_{\mathrm{Cl,charge}} < 0.8\,n_{\mathrm{Cl,ratio}}\\
 0.6\,n_{\mathrm{Cl,charge}} + 0.4\,n_{\mathrm{Cl,ratio}}, & \text{otherwise.}
\end{cases}
```

Report chloride as

```math
\mathrm{Cl\ mg/L} = \max\big(n_{\mathrm{Cl,blend}}\,M_{\mathrm{Cl}}\cdot 1000,\ 0\big).
```

Rationale: charge balance alone can underestimate when inputs miss anions or contain measurement noise; the ratio pathway leverages the stability of major-ion ratios in standard seawater. The adaptive rule prefers the ratio estimate when the two disagree strongly and otherwise uses a conservative blend.

Limitations:

- Assumes composition near standard seawater; large deviations (e.g. hypersaline brines, engineered systems) can bias ratio-derived $n_{\mathrm{Cl}}$.
- Sensitive to analytical errors in high-weight ions (especially Na⁺); propagate uncertainties if available.
- Fixed weights ignore minor temperature/salinity dependence of ionic ratios.
- Fluoride default introduces small uncertainty but typically negligible for overall charge.

### Reference mass per kg

A reference composition in mmol/kg is converted to g/kg using molar masses. Two corrections apply:

1. Replace elemental B by chosen species masses $\mathrm{B(OH)_3},\ \mathrm{B(OH)_4^-}$ per $\alpha_B$.
2. Optionally add a reference alkalinity mass from a chosen `ref_alk_dKH` (default 8.0; 6.2 available for RN compatibility).

Denote the resulting reference total as $\Sigma^{\mathrm{ref}}_{\mathrm{g/kg}}$.

### Relative salinity and fixed-point iteration

Measured mass per liter:

```math
\Sigma_{\mathrm{g/L}} = \sum_j m_j\ \mathrm{[g/L]}.
```

Convert to g/kg via density $ρ$:

```math
\Sigma_{\mathrm{g/kg}} = \frac{\Sigma_{\mathrm{g/L}}}{ρ/1000}.
```

Define relative salinity from the ratio to reference:

```math
SR_{\text{new}} = SR_{\mathrm{REF}}\cdot \frac{\Sigma_{\mathrm{g/kg}}}{\Sigma^{\mathrm{ref}}_{\mathrm{g/kg}}},\quad
SP_{\text{new}}=35\cdot\frac{SR_{\text{new}}}{SR_{\mathrm{REF}}},\quad
SA_{\text{new}}\approx SR_{\text{new}}.
```

Because $ρ$ depends on $SA$ and $CT$, and $SA$ depends on $SP$, iterate to a fixed point:

1. Initialize $SP=35$, $SA = SP\cdot SR_{\mathrm{REF}}/35$.
2. Compute $CT=\mathrm{CT}(SA,T,p)$ and $ρ=\mathrm{ρ}(SA,CT,p)$ via TEOS‑10.
3. Convert all component g/L to g/kg using $ρ$ and sum to $\Sigma_{\mathrm{g/kg}}$.
4. Update $SR$, $SP$, $SA$.
5. Stop when $|SP_{\text{new}}-SP|<\varepsilon$.

Convergence is fast because the density feedback is weak in this range.

### TEOS-10 relations used

- Convert `SP→SA`:

```math
SA \approx SP\cdot\frac{SR_{\mathrm{REF}}}{35}.
```

- Conservative temperature:

```math
CT = \mathrm{CT\_from\_t}(SA, T, p).
```

- Density:

```math
\rho = \mathrm{\rho}(SA, CT, p).
```

- Specific gravity at $t_\mathrm{ref}$:

```math
SG(t_\mathrm{ref}/t_\mathrm{ref})=
\frac{\rho_{\mathrm{sw}}(SA,t_\mathrm{ref},p_\mathrm{ref})}{\rho_{\mathrm{pw}}(SA{=}0,t_\mathrm{ref},p_\mathrm{ref})}.
```

### Component reporting and normalization

For comparability, normalize to `SP = 35`:

```math
  \text{norm\_factor}=\frac{35}{SP},\qquad m^\star=m\cdot\text{norm\_factor}.
```

## References

- Castelao, G. P., and L. Irber (2024): Gibbs Sea Water Oceanographic Toolbox of TEOS‑10 implemented in Rust. Journal of Open Source Software, 9(93), 5988. [https://doi.org/10.21105/joss.05988](https://doi.org/10.21105/joss.05988)
- McDougall, T. J., and P. M. Barker (2011): Getting Started with TEOS‑10 and the Gibbs Seawater (GSW) Oceanographic Toolbox. TEOS‑10 Publication. [https://www.teos-10.org/pubs/Getting_Started.pdf](https://www.teos-10.org/pubs/Getting_Started.pdf)
- IOC, SCOR and IAPSO (2010): The International Thermodynamic Equation of Seawater – 2010 (TEOS‑10) Manual. TEOS‑10 Publication. [https://www.teos-10.org/pubs/TEOS-10_Manual.pdf](https://www.teos-10.org/pubs/TEOS-10_Manual.pdf)
- Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall (2008): The composition of Standard Seawater and the definition of the Reference‑Composition Salinity Scale. Deep‑Sea Research Part I, 55, 50–72. [https://doi.org/10.1016/j.dsr.2008.03.004](https://doi.org/10.1016/j.dsr.2008.03.004)
- UNESCO (1981): The Practical Salinity Scale 1978 (PSS‑78). UNESCO Technical Papers in Marine Science No. 36.

Notes:

- In this tool, SA is taken as approximately equal to SR as per TEOS‑10; the constant SR_REF = 35.16504 g/kg follows Millero et al. (2008).
- The fixed alkalinity split and borate fraction are practical approximations for aquarium contexts rather than strict thermodynamic speciation.

## Build from source

Requires a recent stable Rust toolchain.

- Library only:

```bash
cargo build --release
```

- CLI:

```bash
cargo build --release --features cli
```

Run the CLI:

```bash
target/release/salinity_rs --help
```

## License and acknowledgements

MIT license. See `LICENSE`.

This project builds on TEOS‑10 relations via the excellent `gsw` crate.
