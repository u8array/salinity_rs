#!/usr/bin/env python3
import json
import os
import subprocess
import sys


def main():
    """Minimal example: prints a compact JSON object with required keys.

    This example intentionally omits any error handling — it assumes the
    called binary returns valid JSON containing the required keys.
    """
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cargo_toml = os.path.join(repo_root, "Cargo.toml")

    cmd = [
        "cargo",
        "run",
        "--manifest-path",
        cargo_toml,
        "--quiet",
        "--",
        "--json",
        "--inputs-json",
        json.dumps({
            "na": 11980.0,
            "ca": 357.0,
            "mg": 1246.0,
            "k": 464.0,
            "sr": 6.96,
            "br": 73.2,
            "cl": 19570.0,
            "f": 1.14,
            "s": 814.0,
            "b": 5.57,
            "t_c": 20.0,
            "p_dbar": 0.0,
            "alk_dkh": None,
            "assume_borate": True,
            "default_f_mg_l": 1.296,
            "ref_alk_dkh": 8.0,
            "borate_fraction": None,
            "alk_mg_per_meq": None,
            "return_components": False,
        }),
        "--assumptions-json",
        json.dumps({
            "temp": 20.0,
            "pressure_dbar": 0.0,
            "alkalinity": 8.0,
            "assume_borate": True,
            "default_f_mg_l": 1.296,
            "ref_alk_dkh": 8.0,
            "salinity_norm": 35.0,
            "return_components": False,
            "borate_fraction": None,
            "alk_mg_per_meq": None,
            "rn_compat": False,
        }),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = json.loads(proc.stdout)

    required_keys = [
        "sp",
        "sa_g_per_kg",
        "density_kg_per_m3",
        "sg_20_20",
        "sg_25_25",
    ]

    result = {k: data[k] for k in required_keys}
    print(json.dumps(result, separators=(',', ':')))
    sys.exit(0)


if __name__ == "__main__":
    main()
