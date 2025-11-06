
# salinity_cli_rs

Command-line tool to estimate **Practical Salinity** `SP`, **Absolute Salinity** `SA`, **density** `ρ`, and specific gravities from macro chemical analyses of seawater or reef aquaria.  
Calculations follow TEOS-10 conventions and couple charge/mass balances with density via a simple fixed-point iteration. The library contains hooks for component tables; the current CLI outputs a concise summary (SP, SA, ρ, SG).

## Features

- Compute `SP` and `SA` from elemental/ionic analyses
- TEOS-10 density `ρ(SA,CT,p)` and specific gravity at a reference temperature
- Self-consistent conversion between mg/L and mg/kg via `ρ`
- Chloride estimation from electroneutrality if not measured
- Configurable assumptions: temperature `T`, pressure `p`, alkalinity, borate fraction
- Library support for component tables (mg/L, mg/kg, SP=35 normalization). Note: current CLI prints summary only.

## CLI usage

The binary accepts JSON inputs for measurements and assumptions:

- `--inputs-json <JSON>`: Inline JSON for input values (schema like `Inputs`).
- `--assumptions-json <JSON>`: Optional, adds/overrides assumptions (schema like `Assumptions`).
- `--input <FILE>`: Read a file containing an object with `inputs` and optional `assumptions`. Use `-` for stdin.
- `--json`: Output in JSON format (otherwise human-readable lines for SP/SA/Density/SG).

JSON fields (excerpt):

- Inputs (mg/L unless noted): `na, ca, mg, k, sr, br, s, b, cl` (optional; omit or set `null` for auto-estimate), `f` (optional)
- Conditions: `t_c` (°C, default 20), `p_dbar` (dbar, default 0)
- Alkalinity/options: `alk_dkh` (dKH), `assume_borate` (default true), `borate_fraction`, `ref_alk_dkh` (default 8.0), `alk_mg_per_meq`, `default_f_mg_l` (default 1.296), `return_components` (default false; currently ignored by CLI output)

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

## Install (prebuilt binaries)

Prebuilt binaries for macOS, Linux and Windows are available on the
[Releases](https://github.com/u8array/salinity_cli_rs/releases) page.

### Platform notes

macOS: download the macOS binary for your architecture (Apple Silicon/Intel),
make it executable and put it on your PATH.

```bash
chmod +x ./salinity_teos_10
./salinity_teos_10 --help
```

Linux: download the Linux binary, make it executable and put it on your PATH.

```bash
chmod +x ./salinity_teos_10
./salinity_teos_10 --help
```

Windows: download `salinity_teos_10.exe` and run it from PowerShell or CMD.

```powershell
.\salinity_teos_10.exe --help
```

## Quick start

Example with concentrations in mg/L, `t = 20 °C`, `p = 0 dbar`, and total alkalinity `8 dKH`.

Two invocation methods are supported:

1. Inline JSON (on the command line)

```bash
salinity_teos_10 --inputs-json '{
  "na":10780, "mg":1290, "ca":420, "k":400, "sr":8,
  "br":65, "s":900, "b":4.4,
  "t_c":20, "p_dbar":0, "alk_dkh":8,
  "return_components":true
}'
```

1. File-based (file contains an object with `inputs` and optional `assumptions`)

`inputs.json`

```json
{
  "inputs": {
    "na":10780, "mg":1290, "ca":420, "k":400, "sr":8,
    "br":65, "s":900, "b":4.4,
    "t_c":20, "p_dbar":0, "alk_dkh":8,
    "return_components": true
  }
}
```

Run:

```bash
salinity_teos_10 --input inputs.json
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
2. Compute $CT=\mathrm{CT}(SA,T,p)$, $ρ=\mathrm{ρ}(SA,CT,p)$ via TEOS-10. In the current implementation, $CT$ is approximated by in-situ $t$ (CT ≈ t), which is adequate near 0 dbar.
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
  ext{norm\_factor}=\frac{35}{SP},\qquad m^\star=m\cdot\text{norm\_factor}.
```

## Validation

Compare `SP` against conductivity-derived `SP` and ensure TEOS-10 `ρ(SA,CT,p)` matches expected densities at 20 °C or 25 °C for typical reef mixes.

## References

- McDougall, T. J., and P. M. Barker (2011): Getting Started with TEOS‑10 and the Gibbs Seawater (GSW) Oceanographic Toolbox. TEOS‑10 Publication. [https://www.teos-10.org/pubs/Getting_Started.pdf](https://www.teos-10.org/pubs/Getting_Started.pdf)
- IOC, SCOR and IAPSO (2010): The International Thermodynamic Equation of Seawater – 2010 (TEOS‑10) Manual. TEOS‑10 Publication. [https://www.teos-10.org/pubs/TEOS-10_Manual.pdf](https://www.teos-10.org/pubs/TEOS-10_Manual.pdf)
- Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall (2008): The composition of Standard Seawater and the definition of the Reference‑Composition Salinity Scale. Deep‑Sea Research Part I, 55, 50–72. [https://doi.org/10.1016/j.dsr.2008.03.004](https://doi.org/10.1016/j.dsr.2008.03.004)
- UNESCO (1981): The Practical Salinity Scale 1978 (PSS‑78). UNESCO Technical Papers in Marine Science No. 36.

Notes:

- In this tool, SA is taken as approximately equal to SR as per TEOS‑10; the constant SR_REF = 35.16504 g/kg follows Millero et al. (2008).
- The fixed alkalinity split and borate fraction are practical approximations for aquarium contexts rather than strict thermodynamic speciation.

## Build from source

Requires the latest stable Rust toolchain.

Build locally:

```bash
cargo build --release
```
