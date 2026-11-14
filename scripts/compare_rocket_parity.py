"""Compare rocket stations with native cea==3.3.4; neither solver shares results.

Run in a throwaway environment with numpy, cea==3.3.4 and the installed Rust
wheel. Each case runs in its own process because native failures can abort.
The RP-1311 inputs come from NASA CEA commit 2f79a647, samples/rp1311_examples.inp.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
import warnings
from pathlib import Path

CASES = [
    "8",
    "9",
    "10",
    "11",
    "11-reference",
    "12-eq",
    "12-fr",
    "13",
    "freeze-3",
    "freeze-5",
    "freeze-omit",
    "fac-fr",
    "fac-fr-1",
    "fac-fr-2",
    "fac-fr-3",
    "fac-fr-5",
    "assigned-t",
    "assigned-h",
    "assigned-est",
    "frozen-stop",
    "condensed-stop",
    "freeze-subar",
    "subar-range",
    "pressure-range",
    "fac-sonic",
]
STATE_FIELDS = {
    "temperature_k": ("T", 1),
    "pressure_pa": ("P", 1e5),
    "density_kg_m3": ("density", 1),
    "enthalpy_j_kg": ("enthalpy", 1e3),
    "entropy_j_kg_k": ("entropy", 1e3),
    "gamma_s": ("gamma_s", 1),
    "molecular_weight_kg_kmol": ("M", 1),
}
FLOW_FIELDS = {
    "mach": "Mach",
    "sound_speed_m_s": "sonic_velocity",
    "area_ratio": "ae_at",
    "cstar_m_s": "c_star",
    "thrust_coefficient": "coefficient_of_thrust",
    "isp_m_s": "Isp",
    "isp_vacuum_m_s": "Isp_vacuum",
}


def run_case(case: str) -> dict:
    import cea
    import cea_rs
    import numpy as np

    cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
    assert cea.__version__ == "3.3.4"
    names, weights, temps = ["H2(L)", "O2(L)"], [1, 5.55157], [20.27, 90.17]
    pc, options, only, insert, ions = 53.3172e5, {}, None, [], False
    if case in ("8", "9", "10"):
        options = {"pi_p": [10, 100, 1000], "supar": [25, 50, 75]}
        if case == "8":
            options["subar"] = [1.58]
        else:
            options["iac"] = False
            options["ac_at" if case == "9" else "mdot"] = (
                1.58 if case == "9" else 1333.9
            )
    elif case.startswith("11"):
        names, temps, ions = ["Li(cr)", "F2(L)"], [298.15, 85.02], True
        weights = [6.941, 0.5556 * 37.9968064]
        pc = float(cea.units.psi_to_bar(1000)) * 1e5
        options = {"pi_p": [68.0457], "subar": [10], "supar": [10, 20, 100]}
    elif case.startswith("12"):
        names, weights, temps = ["CH6N2(L)", "N2O4(L)"], [1, 2.5], [298.15] * 2
        pc = float(cea.units.psi_to_bar(1000)) * 1e5
        options = {
            "pi_p": [68.0457],
            "supar": [5, 10, 25, 50, 75, 100, 150, 200],
            "n_frz": 2 if case.endswith("fr") else 0,
        }
        only = [
            "CO",
            "CO2",
            "H",
            "HNO",
            "HNO2",
            "HO2",
            "H2",
            "H2O",
            "H2O2",
            "N",
            "NO",
            "NO2",
            "N2",
            "N2O",
            "O",
            "OH",
            "O2",
            "HCO",
            "NH",
            "CH4",
            "NH2",
            "NH3",
            "H2O(L)",
            "C(gr)",
        ]
    elif case == "13":
        names = ["N2H4(L)", "Be(a)", "H2O2(L)"]
        weights, temps = [0.67 * 0.8, 0.67 * 0.2, 0.33], [298.15] * 3
        pc = float(cea.units.psi_to_bar(3000)) * 1e5
        options, insert = {"pi_p": [3, 10, 30, 300]}, ["BeO(L)"]
    elif case == "freeze-3":
        options = {"pi_p": [10, 100, 1000], "supar": [25, 50, 75], "n_frz": 3}
    elif case == "freeze-5":
        options = {"supar": [2, 5, 10, 25, 50], "n_frz": 5}
    elif case == "freeze-omit":
        options = {
            "pi_p": [10, 2, 100],
            "subar": [2],
            "supar": [1.5, 25, 10, 50],
            "n_frz": 3,
        }
    elif case.startswith("fac-fr"):
        options = {
            "iac": False,
            "ac_at": 1.58,
            "supar": [5, 25, 50],
            "n_frz": 4 if case == "fac-fr" else int(case.rsplit("-", 1)[1]),
        }
    elif case == "assigned-t":
        options = {"tc": 3200, "pi_p": [10], "subar": [1.5], "supar": [5, 25]}
    elif case == "assigned-h":
        options = {"hc": 1e6, "pi_p": [10], "supar": [5, 25]}
    elif case == "assigned-est":
        options = {"tc_est": 3200, "supar": [5, 25]}
    elif case == "frozen-stop":
        names, weights, temps = ["N2"], [1], [300]
        options = {"tc": 500, "pi_p": [2, 5, 10, 100], "n_frz": 1}

    elif case == "condensed-stop":
        names, weights, temps = ["CO"], [1], [298.15]
        options = {"tc": 900, "pi_p": [2, 3, 5, 10, 100], "n_frz": 1}
    elif case == "freeze-subar":
        options = {"pi_p": [10], "subar": [1.58], "supar": [5, 25], "n_frz": 1}
    elif case == "subar-range":
        options = {"subar": [1.001, 1.01, 2, 10, 100]}
    elif case == "pressure-range":
        options = {"pi_p": [1.01, 1.1, 2, 10, 100, 1000]}

    elif case == "fac-sonic":
        options = {"iac": False, "ac_at": 1.0, "supar": [5, 25]}

    reac = cea.Mixture(names, ions=ions)
    prod = (
        cea.Mixture(only, ions=ions)
        if only
        else cea.Mixture(names, products_from_reactants=True, ions=ions)
    )
    native_solver = cea.RocketSolver(prod, reactants=reac, ions=ions, insert=insert)
    native = cea.RocketSolution(native_solver)
    native_options = dict(options)
    if "hc" in native_options:
        native_options["hc"] /= cea.R
    elif "tc" not in native_options:
        native_options["hc"] = (
            reac.calc_property(cea.ENTHALPY, np.array(weights), np.array(temps)) / cea.R
        )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        start = time.perf_counter()
        native_solver.solve(
            native, np.array(weights, dtype=float), pc / 1e5, **native_options
        )
        native_time = time.perf_counter() - start
    result = {
        "case": case,
        "inputs": {
            "names": names,
            "weights": weights,
            "temperatures": temps,
            "pressure_pa": pc,
            "only": only,
            "insert": insert,
            "ions": ions,
            "options": options,
        },
        "native_converged": bool(native.converged),
        "native_error": int(native.last_error),
        "native_warnings": [str(w.message) for w in caught],
        "native_seconds": native_time,
    }
    # Keep every native station as a reusable golden, even if Rust cannot solve.
    native_rows = []
    for i in range(native.num_pts):
        row = {
            key: float(getattr(native, prop)[i]) * scale
            for key, (prop, scale) in STATE_FIELDS.items()
        }
        row.update(
            {key: float(getattr(native, prop)[i]) for key, prop in FLOW_FIELDS.items()}
        )
        row["velocity_m_s"] = float(native.Mach[i] * native.sonic_velocity[i])
        native_rows.append(row)
    native_finite = all(math.isfinite(v) for row in native_rows for v in row.values())
    result["native_finite"] = native_finite
    result["native"] = [
        {k: v if math.isfinite(v) else None for k, v in row.items()}
        for row in native_rows
    ]
    try:
        rust_names = names
        if case == "11-reference":
            # Li(cr)'s NASA record declares a reversed first interval and is omitted
            # from this crate. Its standard-state reference enthalpy is exactly zero.
            rust_names = [
                cea_rs.Reactant(
                    "Li-reference",
                    formula={"Li": 1},
                    molecular_weight=6.941,
                    enthalpy_j_mol=0,
                ),
                "F2(L)",
            ]
            result["rust_reactant_override"] = (
                "Li(cr) replaced by formula Li=1, M=6.941, H=0 J/mol at 298.15 K"
            )
        solver = cea_rs.Solver(rust_names, only=only, insert=insert, ions=ions)
        start = time.perf_counter()
        rust = solver.rocket(weights, temps, pc, **options)
        result["rust_seconds"] = time.perf_counter() - start
        result["rust_complete"] = rust.complete
        result["rust_stop_reason"] = rust.stop_reason
        rows, differences = [], []
        for i, station in enumerate(rust.stations):
            row = {key: float(getattr(station.state, key)) for key in STATE_FIELDS}
            row.update({key: float(getattr(station, key) or 0) for key in FLOW_FIELDS})
            row["velocity_m_s"] = station.velocity_m_s
            rows.append(row)
            if native_finite and i < len(native_rows):
                differences.append(
                    {
                        k: abs(v - native_rows[i][k])
                        / max(abs(native_rows[i][k]), 1e-30)
                        if native_rows[i][k]
                        else abs(v)
                        for k, v in row.items()
                    }
                )
        result.update(
            rust=rows,
            relative_differences=differences,
            worst_relative_difference=(
                max((v for d in differences for v in d.values()), default=0)
                if native_finite
                else None
            ),
        )
        result["passed"] = (
            native_finite
            and len(rows) == len(native_rows)
            and result["worst_relative_difference"] <= 1e-4
            and (native.converged == rust.complete)
        )
    except (ValueError, RuntimeError) as exc:
        result.update(rust_error=f"{type(exc).__name__}: {exc}", passed=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=CASES)
    args = parser.parse_args()
    if args.case:
        report = run_case(args.case)
    else:
        reports = []
        for case in CASES:
            output = args.output.with_name(args.output.stem + "-" + case + ".json")
            proc = subprocess.run(
                [sys.executable, __file__, "--case", case, "--output", str(output)],
                check=False,
                capture_output=True,
                text=True,
                timeout=180,
            )
            if proc.returncode:
                reports.append(
                    {
                        "case": case,
                        "passed": False,
                        "process_error": proc.stdout + proc.stderr,
                    }
                )
            else:
                reports.append(json.loads(output.read_text()))
            r = reports[-1]
            print(
                case,
                r["passed"],
                r.get(
                    "worst_relative_difference",
                    r.get("rust_error", r.get("process_error")),
                ),
                flush=True,
            )
        report = {
            "oracle": "cea==3.3.4",
            "nasa_commit": "2f79a647e85737742fb08a127da00d1ea1d069da",
            "tolerance": 1e-4,
            "scale": "absolute relative difference; absolute error at zero",
            "cases": reports,
            "all_passed": all(r["passed"] for r in reports),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if not args.case and not report["all_passed"]:
        raise SystemExit("Rocket comparisons include failures; see the recorded cases")


if __name__ == "__main__":
    main()
