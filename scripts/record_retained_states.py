"""Record NASA 3.3.4 failures, false positives and retained fields.

Requires numpy, cea==3.3.4, and an installed cea_rs wheel for the shared test
helpers. Expected values are read only from the independently installed cea.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import runpy
import sys
import warnings

import numpy as np


def clean(value):
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, list):
        return [clean(item) for item in value]
    if isinstance(value, dict):
        return {key: clean(item) for key, item in value.items()}
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installed", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("tests/data/nasa/retained-3.3.4.json"))
    args = parser.parse_args()
    sys.argv = [__file__, "--installed", str(args.installed.resolve())]
    helpers = runpy.run_path(str(Path(__file__).resolve().parents[1] / "tests/test_nasa_api.py"), run_name="record_retained")
    api = helpers["native"]
    if api is None:
        raise SystemExit("Recording requires cea==3.3.4")
    api.set_log_level(api.LOG_NONE)
    inputs = {
        "hg-hp": ("condensed", ["Hg", "HP", -15551.803948663037 / api.R, [1.0, 20.0]], {}),
        "water-transition": ("condensed", ["H2O", "TP", 273.15, [1.0, 1.0]], {}),
        "nacn-transition": ("condensed", ["NaCN", "TP", 288.5, [1.0, 1.0]], {}),
        "nacn-300": ("condensed", ["NaCN", "TP", 300.0, [1.0, 1.0]], {}),
        "br2-supersaturation": ("condensed", ["Br2", "TP", 262.0, [1.0, 20.0]], {}),
        "air-jet-a": ("air_jet_a", [], {}),
        "ap-al-nan": ("ap_al", [], {}),
        "ap-al-retained-transport": ("ap_al", [], {"warm_start": True}),
        "water-reuse": ("reused_failure", [], {}),
        "argon-low-temperature": ("argon_below_range", [], {}),
        "fac-early-freeze": ("rp_rocket", [9, 3], {}),
        "shock-last-valid": ("shock", [], {"speed": 1000.0}),
        "shock-retained-reflection": ("shock", [], {"speed": 1000.0, "reflected": True}),
        "nacn-reduced": ("condensed", ["NaCN", "TP", 300.0, [1.0, 1.0]], {"pressure": 1e20}),
        "shock-reflected-failure": ("shock_edge", ["reflected-failure"], {}),
        "shock-incident-failure": ("shock_edge", ["incident-failure"], {}),
        "shock-atomic-transport": ("shock_edge", ["atomic-transport"], {}),
        "shock-incident-cap": ("shock_edge", ["incident-cap"], {}),
        "shock-reflected-cap": ("shock_edge", ["reflected-cap"], {}),
        "shock-singular-last-valid": ("shock_edge", ["last-valid"], {}),
        "shock-ion-last-valid": ("shock_edge", ["ion-last-valid"], {}),
        "detonation-nonfinite": ("nonfinite_detonation", [], {}),
        "detonation-hot-nonfinite": ("nonfinite_detonation", [], {"hot": True}),
        "initial-guess-failure": ("initial_guess_failure", [], {}),
    }
    cases = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for key, (name, positional, keyword) in inputs.items():
            solution, _, _ = helpers[name](api, *positional, **keyword)
            if key == "shock-ion-last-valid":
                # Compare status/validity: singular ion iterates differ numerically across solvers.
                for field in ("T", "P", "velocity", "sonic_velocity"):
                    value = getattr(solution, field)
                    assert np.all(np.isfinite(value)) and np.all(value > 0.0)
                cases[key] = {
                    "converged": solution.converged,
                    "last_error": solution.last_error,
                    "values": {},
                    "diagnostic_values": {
                        field: clean(getattr(solution, field))
                        for field in ("T", "P", "P21", "T21")
                    },
                }
                continue
            fields = helpers["THERMO"] + helpers["TRANSPORT"]
            if name in ("rp_rocket", "air_jet_a"):
                fields += helpers["PERFORMANCE"] + ("n", "nj")
            elif name in ("condensed", "ap_al", "reused_failure", "argon_below_range", "initial_guess_failure"):
                fields += ("volume", "n", "nj", "ln_nj")
            elif name == "nonfinite_detonation":
                fields += ("Mach", "velocity", "sonic_velocity", "P_P1", "T_T1", "M_M1", "rho_rho1", "n", "nj", "ln_nj")
            else:
                fields += helpers["SHOCK_FIELDS"]
            values = {field: clean(getattr(solution, field)) for field in fields}
            values.update({field: clean(getattr(solution, field)) for field in ("mass_fractions", "mole_fractions")})
            cases[key] = {"converged": solution.converged, "last_error": solution.last_error, "values": values}
    report = {
        "native_version": api.__version__,
        "nasa_commit": "2f79a647e85737742fb08a127da00d1ea1d069da",
        "inputs": inputs,
        "cases": cases,
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
