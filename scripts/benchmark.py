"""Measure Rust release-mode and native CEA with the same physical inputs.

This measures single-process elapsed solve time, not operating-system disk-cache
cold starts. The native solver is not warm-started from its previous solution:
a fresh solution object is allocated per call, as in Hybrid's current consumer.
Construction/first-call times are separate from repeated warm solves.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--installed", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--repeats", type=int, default=50)
parser.add_argument("--cold-backend", choices=["native", "rust"])
parser.add_argument("--cold-mode", default="equilibrium")
args = parser.parse_args()
sys.path.insert(0, str(args.installed.resolve()))
if args.cold_backend:
    start = time.perf_counter_ns()
    if args.cold_backend == "native":
        os.environ["CEA_SKIP_INIT"] = "1"
        import cea
        import numpy as np

        cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
    else:
        import cea_rs
    imported = time.perf_counter_ns()
    weights = [1.0, 6.0]
    if args.cold_backend == "native":
        reac = cea.Mixture(["H2", "O2"])
        prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
        solver = (
            cea.EqSolver(prod, reactants=reac)
            if args.cold_mode == "TP"
            else cea.RocketSolver(prod, reactants=reac, transport=False)
        )
        hc = (
            reac.calc_property(cea.ENTHALPY, np.array(weights), np.array([298.15] * 2))
            / cea.R
        )
        sol = (
            cea.EqSolution(solver)
            if args.cold_mode == "TP"
            else cea.RocketSolution(solver)
        )
    else:
        solver = cea_rs.Solver(["H2", "O2"])
    constructed = time.perf_counter_ns()
    if args.cold_backend == "native":
        if args.cold_mode == "TP":
            solver.solve(sol, cea.TP, 3000.0, 30.0, np.array(weights))
        else:
            solver.solve(
                sol,
                np.array(weights),
                30.0,
                supar=10.0,
                hc=hc,
                iac=True,
                n_frz={"equilibrium": 0, "frozen_at_throat": 2}[args.cold_mode],
            )
        if not sol.converged or sol.last_error:
            raise RuntimeError("Native process-cold solve failed")
    elif args.cold_mode == "TP":
        solver.solve("TP", 3000.0, 3e6, weights)
    else:
        solver.rocket(weights, [298.15] * 2, 3e6, 10.0, flow=args.cold_mode)
    finished = time.perf_counter_ns()
    print(
        json.dumps(
            {
                "import_init_ms": (imported - start) / 1e6,
                "construction_ms": (constructed - imported) / 1e6,
                "first_solve_ms": (finished - constructed) / 1e6,
            }
        )
    )
    raise SystemExit(0)

import cea
import cea_rs
import numpy as np

cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)


def measure(call: Callable[[], object], repeats: int) -> dict[str, float]:
    """Measure identical calls after five discarded warmup solves."""
    for _ in range(5):
        call()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        call()
        samples.append((time.perf_counter_ns() - start) / 1e6)
    return {
        "median_ms": statistics.median(samples),
        "min_ms": min(samples),
        "max_ms": max(samples),
    }


rows = []
for names, weights, temps in [
    (["H2", "O2"], [1.0, 6.0], [298.15] * 2),
    (["O2(L)", "RP-1"], [2.6, 1.0], [90.17, 298.15]),
    (["N2O4(L)", "N2H4(L)", "C2H8N2(L),UDMH"], [1.6, 0.5, 0.5], [298.15] * 3),
]:
    start = time.perf_counter_ns()
    r = cea.Mixture(names)
    p = cea.Mixture(names, products_from_reactants=True)
    native_eq = cea.EqSolver(p, reactants=r)
    native_rocket = cea.RocketSolver(p, reactants=r, transport=False)
    hc = r.calc_property(cea.ENTHALPY, np.array(weights), np.array(temps)) / cea.R
    native_construct = (time.perf_counter_ns() - start) / 1e6
    start = time.perf_counter_ns()
    rust = cea_rs.Solver(names)
    rust_construct = (time.perf_counter_ns() - start) / 1e6
    for mode in ["TP", "equilibrium", "frozen_at_chamber", "frozen_at_throat"]:

        def native_call(
            mode: str = mode,
            native_eq: cea.EqSolver = native_eq,
            native_rocket: cea.RocketSolver = native_rocket,
            weights: list[float] = weights,
            hc: float = hc,
        ) -> None:
            """Allocate and verify one native CEA solution."""
            if mode == "TP":
                sol = cea.EqSolution(native_eq)
                native_eq.solve(sol, cea.TP, 3000.0, 30.0, np.array(weights))
            else:
                sol = cea.RocketSolution(native_rocket)
                native_rocket.solve(
                    sol,
                    np.array(weights),
                    30.0,
                    supar=10.0,
                    hc=hc,
                    iac=True,
                    n_frz={
                        "equilibrium": 0,
                        "frozen_at_chamber": 1,
                        "frozen_at_throat": 2,
                    }[mode],
                )
            if not sol.converged or sol.last_error:
                raise RuntimeError("Native benchmark solve failed")

        def rust_call(
            mode: str = mode,
            rust: cea_rs.Solver = rust,
            weights: list[float] = weights,
            temps: list[float] = temps,
        ) -> cea_rs.State | cea_rs.RocketResult:
            """Run the matching validated Rust solve."""
            return (
                rust.solve("TP", 3000.0, 3e6, weights)
                if mode == "TP"
                else rust.rocket(weights, temps, 3e6, 10.0, flow=mode)
            )

        first = {}
        for name, call in [("native", native_call), ("rust", rust_call)]:
            start = time.perf_counter_ns()
            call()
            first[name] = (time.perf_counter_ns() - start) / 1e6
        # Alternate which implementation is measured first across problem types.
        order = [("native", native_call), ("rust", rust_call)]
        if len(rows) % 2:
            order.reverse()
        measurements = {name: measure(call, args.repeats) for name, call in order}
        rows.append(
            {
                "reactants": names,
                "mode": mode,
                "construction_ms": {"native": native_construct, "rust": rust_construct},
                "first_call_ms": first,
                "warm": measurements,
                "rust_over_native_warm_median": measurements["rust"]["median_ms"]
                / measurements["native"]["median_ms"],
            }
        )
        print(
            names,
            mode,
            "rust/native",
            rows[-1]["rust_over_native_warm_median"],
            flush=True,
        )
cold_rows = []
for mode in ["TP", "equilibrium", "frozen_at_throat"]:
    for backend in ["native", "rust"]:
        samples = []
        for _ in range(5):
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--installed",
                    str(args.installed.resolve()),
                    "--output",
                    str(args.output),
                    "--cold-backend",
                    backend,
                    "--cold-mode",
                    mode,
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            samples.append(json.loads(completed.stdout))
        cold_rows.append({"mode": mode, "backend": backend, "samples_ms": samples})
report = {
    "python": platform.python_version(),
    "platform": platform.platform(),
    "native": cea.__version__,
    "rust": cea_rs.__version__,
    "profile": "cargo --release",
    "warmups": 5,
    "repeats": args.repeats,
    "cache_conditions": "Same process and installed native libraries. First-call and construction recorded; warm repeats after five discarded calls. OS filesystem caches not flushed. Fresh native solution per call; no prior-solution warm start. Full validation remains enabled in Rust.",
    "process_cold_cases": cold_rows,
    "process_cold_conditions": "Five fresh Python processes per backend/mode, H2/O2 1:6, 3 MPa, Ae/At=10; import/init, construction and first solve timed separately. OS filesystem caches remain warm/uncontrolled; no disk-cold claim.",
    "cases": rows,
}
args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
