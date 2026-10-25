"""Compare an installed cea_rs wheel with independently installed NASA CEA.

Run with a Python interpreter containing cea==3.3.4 and numpy. Pass --installed
pointing to a fresh wheel installation outside the source tree. Both bindings
are independent native libraries; no oracle output enters the Rust solve.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--installed", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--section", required=True)
args = parser.parse_args()
sys.path.insert(0, str(args.installed.resolve()))
import cea
import cea_rs
import numpy as np

cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
if cea.__version__ != "3.3.4":
    raise RuntimeError(
        f"Require the recorded native oracle cea==3.3.4, got {cea.__version__}"
    )

FIELDS = {
    "temperature_k": ("T", 1.0),
    "pressure_pa": ("P", 1e5),
    "density_kg_m3": ("density", 1.0),
    "enthalpy_j_kg": ("enthalpy", 1000.0),
    "energy_j_kg": ("energy", 1000.0),
    "entropy_j_kg_k": ("entropy", 1000.0),
    "cp_equilibrium_j_kg_k": ("cp_eq", 1000.0),
    "cp_frozen_j_kg_k": ("cp_fr", 1000.0),
    "gamma_s": ("gamma_s", 1.0),
    "molecular_weight_kg_kmol": ("M", 1.0),
}
rows: list[dict[str, object]] = []


def compare(
    label: str,
    native: dict[str, float],
    rust: dict[str, float],
    tolerance: float = 3e-4,
) -> None:
    """Record every field and fail on nonfinite or excessive scaled difference."""
    delta = {
        k: abs(rust.get(k, 0.0) - native.get(k, 0.0))
        / max(abs(native.get(k, 0.0)), 1.0)
        for k in native.keys() | rust.keys()
    }
    passed = all(math.isfinite(v) and v <= tolerance for v in delta.values())
    rows.append(
        {
            "case": label,
            "native": native,
            "rust": rust,
            "scaled_differences": delta,
            "tolerance": tolerance,
            "passed": passed,
        }
    )
    print(
        label,
        "PASS" if passed else "FAIL",
        f"max={max(delta.values()):.6g}",
        flush=True,
    )


def mixtures(
    names: list[str | dict[str, object]],
    only: list[str] | None = None,
    ions: bool = False,
) -> tuple[cea.Mixture, cea.Mixture, cea_rs.Solver]:
    """Build matching declared reactant/product models in each implementation."""
    native_names = [
        cea.Reactant(
            name=r["name"],
            formula=r["formula"],
            molecular_weight=r["molecular_weight"],
            enthalpy=r["enthalpy_j_mol"],
            enthalpy_units="J/mol",
            temperature=298.15,
        )
        if isinstance(r, dict)
        else r
        for r in names
    ]
    rust_names = [cea_rs.Reactant(**r) if isinstance(r, dict) else r for r in names]
    reactants = cea.Mixture(native_names, ions=ions)
    products = (
        cea.Mixture(only, ions=ions)
        if only
        else cea.Mixture(native_names, products_from_reactants=True, ions=ions)
    )
    return reactants, products, cea_rs.Solver(rust_names, only=only, ions=ions)


# Multiple elements, pressures and temperatures, with unconstrained product sets.
for case_index, (names, weights) in enumerate(
    [(["H2", "O2"], [1.0, 8.0]), (["CH4", "O2", "N2"], [1.0, 4.0, 13.0])]
):
    if args.section != f"equilibrium-{case_index}":
        continue
    reac, prod, rust = mixtures(names)
    native_solver = cea.EqSolver(prod, reactants=reac)
    for t, p in [(1000.0, 1e5), (2500.0, 1e6), (3500.0, 5e6)]:
        native = cea.EqSolution(native_solver)
        native_solver.solve(native, cea.TP, t, p / 1e5, np.array(weights))
        if not native.converged or native.last_error:
            raise RuntimeError("Native oracle did not converge")
        result = rust.solve("TP", t, p, weights)
        compare(
            f"TP:{names}:{t}:{p}",
            {
                k: float(getattr(native, field)) * scale
                for k, (field, scale) in FIELDS.items()
            },
            {k: getattr(result, k) for k in FIELDS},
        )
        tp_result = result
        # Compare all six modes against the oracle, not solely Rust round trips.
        if t == 2500.0:
            state = {
                "HP": native.enthalpy * 1000,
                "SP": native.entropy * 1000,
                "TV": t,
                "UV": native.energy * 1000,
                "SV": native.entropy * 1000,
            }
            for mode, first in state.items():
                second = p if mode.endswith("P") else native.volume
                oracle = cea.EqSolution(native_solver)
                native_first = first if mode.startswith("T") else first / cea.R
                native_solver.solve(
                    oracle,
                    getattr(cea, mode),
                    native_first,
                    second / 1e5 if mode.endswith("P") else second,
                    np.array(weights),
                )
                if not oracle.converged or oracle.last_error:
                    raise RuntimeError(f"Native {mode} oracle did not converge")
                result = rust.solve(mode, first, second, weights)
                compare(
                    f"{mode}:{names}",
                    {
                        k: float(getattr(oracle, field)) * scale
                        for k, (field, scale) in FIELDS.items()
                    },
                    {k: getattr(result, k) for k in FIELDS},
                )
        compare(
            f"composition:{names}:{t}:{p}",
            {k: float(v) for k, v in native.mole_fractions.items()},
            tp_result.mole_fractions,
            3e-5,
        )

# Condensed and ions: explicitly declared pools, identical in both solvers.
for names, weights, t, only, ions in [
    (["CO"], [1.0], 800.0, ["CO", "CO2", "C", "O", "O2", "C(gr)"], False),
    (["H2"], [1.0], 10000.0, ["H", "H2", "H+", "e-"], True),
]:
    if args.section != ("ions" if ions else "condensed"):
        continue
    reac, prod, rust = mixtures(names, only, ions)
    solver = cea.EqSolver(prod, reactants=reac, ions=ions)
    native = cea.EqSolution(solver)
    solver.solve(native, cea.TP, t, 1.0, np.array(weights))
    if not native.converged or native.last_error:
        raise RuntimeError("Native condensed/ion oracle did not converge")
    result = rust.solve("TP", t, 1e5, weights)
    compare(
        f"special:{names}:{ions}",
        {
            k: float(getattr(native, field)) * scale
            for k, (field, scale) in FIELDS.items()
        },
        {k: getattr(result, k) for k in FIELDS},
    )
    compare(
        f"special-composition:{names}:{ions}",
        {k: float(v) for k, v in native.mole_fractions.items()},
        result.mole_fractions,
        3e-5,
    )

# Infinite-area rocket performance, including native frozen-location controls.
ABS = [
    {
        "name": "Acrylonitrile",
        "formula": {"C": 3, "H": 3, "N": 1},
        "molecular_weight": 53.0626,
        "enthalpy_j_mol": 140000.0,
    },
    {
        "name": "Butadiene",
        "formula": {"C": 4, "H": 6},
        "molecular_weight": 54.0904,
        "enthalpy_j_mol": 90500.0,
    },
    {
        "name": "Styrene",
        "formula": {"C": 8, "H": 8},
        "molecular_weight": 104.1491,
        "enthalpy_j_mol": 103400.0,
    },
]
for case_index, (names, weights, temperatures) in enumerate(
    [
        (["H2", "O2"], [1.0, 6.0], [298.15, 298.15]),
        (["O2(L)", "RP-1"], [2.6, 1.0], [90.17, 298.15]),
        (["N2O", "Paraffin"], [6.0, 1.0], [300.0, 298.15]),
        (["N2O", *ABS], [6.0, 0.3992, 0.4737, 0.1271], [300.0, 298.15, 298.15, 298.15]),
        (["N2O4(L)", "N2H4(L)", "C2H8N2(L),UDMH"], [1.6, 0.5, 0.5], [298.15] * 3),
    ]
):
    if args.section != f"rocket-{case_index}":
        continue
    reac, prod, rust = mixtures(names)
    native_solver = cea.RocketSolver(prod, reactants=reac, transport=False)
    hc = (
        reac.calc_property(cea.ENTHALPY, np.array(weights), np.array(temperatures))
        / cea.R
    )
    for pc, area in [(1e6, 3.68), (5e6, 20.0)]:
        for flow, freeze in [
            ("equilibrium", 0),
            ("frozen_at_chamber", 1),
            ("frozen_at_throat", 2),
        ]:
            native = cea.RocketSolution(native_solver)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                start = time.perf_counter()
                native_solver.solve(
                    native,
                    np.array(weights),
                    pc / 1e5,
                    supar=area,
                    hc=hc,
                    iac=True,
                    n_frz=freeze,
                )
                native_seconds = time.perf_counter() - start
            if not native.converged or native.last_error or caught:
                raise RuntimeError(
                    f"Native rocket oracle failed: {names} {flow}: {[str(w.message) for w in caught]}"
                )
            start = time.perf_counter()
            result = rust.rocket(weights, temperatures, pc, area, flow=flow)
            rust_seconds = time.perf_counter() - start
            for i, station in enumerate([result.chamber, result.throat, result.exit]):
                n = {
                    k: float(getattr(native, field)[i]) * scale
                    for k, (field, scale) in FIELDS.items()
                }
                r = {k: getattr(station.state, k) for k in FIELDS}
                # Frozen station cp_eq in NASA is retained equilibrium diagnostic;
                # compare actual frozen Cp, gamma, velocity and state independently.
                if flow != "equilibrium":
                    n.pop("cp_equilibrium_j_kg_k")
                    r.pop("cp_equilibrium_j_kg_k")
                n["velocity_m_s"] = float(native.Mach[i] * native.sonic_velocity[i])
                r["velocity_m_s"] = station.velocity_m_s
                compare(f"rocket:{names}:{pc}:{area}:{flow}:{i}", n, r, 1e-4)
            compare(
                f"rocket-performance:{names}:{pc}:{area}:{flow}",
                {
                    "cstar_m_s": float(native.c_star[0]),
                    "vacuum_velocity_m_s": float(native.Isp_vacuum[-1]),
                },
                {
                    "cstar_m_s": result.cstar_m_s,
                    "vacuum_velocity_m_s": result.vacuum_velocity_m_s,
                },
                1e-4,
            )
            rows[-1]["timing_single_call_seconds"] = {
                "native": native_seconds,
                "rust": rust_seconds,
            }

loaded_native = Path(importlib.import_module("cea.lib.libcea").__file__).resolve()
if not loaded_native.is_file() or loaded_native.suffix not in (".so", ".pyd", ".dylib"):
    raise RuntimeError("Cannot identify loaded native CEA extension")
files = [loaded_native, *list((Path(cea.__file__).parent / "data").glob("*.lib"))]
if len(files) < 2:
    raise RuntimeError("Cannot identify packaged native CEA data")
rust_data = (
    Path(__file__).resolve().parents[1] / "vendor/thermo-species/data/thermo.inp"
)
rust_data_hash = hashlib.sha256(rust_data.read_bytes()).hexdigest()
if rust_data_hash != cea_rs.NASA_THERMO_SHA256:
    raise RuntimeError("Installed Rust wheel and source NASA data identity differ")

report = {
    "native_version": cea.__version__,
    "native_library_version": cea.lib_version(),
    "rust_version": cea_rs.__version__,
    "python": platform.python_version(),
    "platform": platform.platform(),
    "rust_module": cea_rs.__file__,
    "rust_data_sha256": rust_data_hash,
    "native_loaded_extension": str(loaded_native),
    "native_files_sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files
    },
    "scale_definition": "abs(rust-native)/max(abs(native), 1 SI unit)",
    "cases": rows,
    "all_passed": all(row["passed"] for row in rows),
    "scope": "Finite cases only; not full CEA parity, flight qualification, or a performance benchmark.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
if not report["all_passed"]:
    raise SystemExit("Native comparison tolerance exceeded; see JSON evidence")
