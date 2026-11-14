"""Record native cea==3.3.4 equilibrium regressions; never execute reference scripts.

Run in a throwaway project containing cea==3.3.4 and numpy:
    python scripts/record_equilibrium_parity.py --nasa-source /tmp/nasa-cea
The reference checkout must be commit 2f79a647e85737742fb08a127da00d1ea1d069da.
Only literal sample inputs are read from it. Output is consumed by tests/solver.rs.
"""

import argparse
import ast
from pathlib import Path

import cea
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--nasa-source", type=Path, required=True)
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "tests/fixtures/equilibrium-native.txt")
args = parser.parse_args()
if cea.__version__ != "3.3.4":
    raise RuntimeError("The recorded oracle must be cea==3.3.4")
cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
FIELDS = [("T", 1), ("P", 1e5), ("density", 1), ("enthalpy", 1000),
          ("entropy", 1000), ("cp_eq", 1000), ("gamma_s", 1), ("M", 1)]
rows = []


def literal(example, name):
    source = args.nasa_source / f"source/bind/python/cea/samples/rp1311/example{example}.py"
    for node in ast.parse(source.read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(f"Missing literal {name} in {source}")


def record(label, names, weights, mode, first, second, solution, *, only=(), omit=(), trace=0, temps=()):
    if not solution.converged or solution.last_error:
        raise RuntimeError(f"Native oracle failed: {label}, error={solution.last_error}")
    numbers = lambda values: ";".join(str(float(v)) for v in values)
    fractions = ";".join(f"{name}={float(x)}" for name, x in solution.mole_fractions.items() if x > 1e-12)
    rows.append("|".join([label, ";".join(names), numbers(weights), mode, str(first), str(second),
                          str(trace), ";".join(only), ";".join(omit), numbers(temps),
                          numbers(float(getattr(solution, f)) * scale for f, scale in FIELDS), fractions]))


# HP chambers and independent SP solves at native rocket exit pressures/entropies.
for names, weights, pressures, temps in [
    (["NH4CLO4(I)", "AL(cr)"], [4, 1], [.7, 3, 7, 20], [300., 300.]),
    (["NH4CLO4(I)", "AL(cr)"], [6, 1], [.7, 3, 7, 20], [300., 300.]),
    (["NH4CLO4(I)", "AL(cr)", "C4H6,butadiene"], [68, 18, 14], [.7, 3, 7, 20], [300.] * 3),
    (["Air", "Jet-A(L)"], [15, 1], [7, 20], [298.15] * 2),
    (["Air", "Jet-A(L)"], [15, 1], [7, 20], [300.] * 2),
    (["O2(L)", "H2(L)"], [3, 1], [3], [90.17, 20.27]),
]:
    reactants = cea.Mixture(names)
    products = cea.Mixture(names, products_from_reactants=True)
    solver = cea.EqSolver(products, reactants=reactants)
    rocket_solver = cea.RocketSolver(products, reactants=reactants)
    weights = np.array(weights, dtype=float)
    h = reactants.calc_property(cea.ENTHALPY, weights, np.array(temps))
    for pressure in pressures:
        label = f"{names}:{weights.tolist()}:{temps}K:{pressure}MPa"
        solution = cea.EqSolution(solver)
        solver.solve(solution, cea.HP, h / cea.R, pressure * 10, weights)
        record(label, names, weights, "HP", h, pressure * 1e6, solution, temps=temps)
        for area in [15., 50., 200.]:
            rocket = cea.RocketSolution(rocket_solver)
            rocket_solver.solve(rocket, weights, pressure * 10, supar=area, hc=h / cea.R)
            if not rocket.converged or rocket.last_error:
                raise RuntimeError(f"Native rocket failed: {label}, area={area}")
            # Re-solve independently to verify Solver::solve rather than a nozzle warm start.
            solution = cea.EqSolution(solver)
            entropy = float(rocket.entropy[0]) * 1000
            p = float(rocket.P[-1]) * 1e5
            solver.solve(solution, cea.SP, entropy / cea.R, p / 1e5, weights)
            record(f"{label}:area={area}", names, weights, "SP", entropy, p, solution)

# RP-1311 input literals are data, not imported/executed Python.
for example in [1, 2, 3, 4, 5, 14]:
    only, omit, trace = [], [], 0
    if example in (1, 2):
        names = ["H2", "Air"]
        only = literal(example, "prod_names")
    elif example in (3, 4):
        names = ["Air", "C7H8(L)", "C8H18(L),n-octa"]
        omit, trace = literal(3, "omit_names"), 1e-15
    elif example == 5:
        names = ["NH4CLO4(I)", "CHOS-Binder", "AL(cr)", "MgO(cr)", "H2O(L)"]
        omit = literal(5, "omit_names")
    else:
        names = ["H2(L)", "O2(L)"]
    native_names = [cea.Reactant(name=n, formula={"C": 1., "H": 1.86955, "O": .031256, "S": .008415},
                               molecular_weight=14.6652984484, enthalpy=-2999.082 * 4.184,
                               enthalpy_units="J/mol", temperature=298.15) if n == "CHOS-Binder" else n for n in names]
    reactants = cea.Mixture(native_names)
    products = cea.Mixture(only) if only else cea.Mixture(native_names, products_from_reactants=True, omit=omit)
    solver = cea.EqSolver(products, reactants=reactants, trace=trace)
    inputs = []
    if example in (1, 2):
        fuel = reactants.moles_to_weights(np.array([1., 0.]))
        oxidant = reactants.moles_to_weights(np.array([0., 1.]))
        if example == 1:
            for ratio in [1., 1.5]:
                of = reactants.chem_eq_ratio_to_of_ratio(oxidant, fuel, ratio)
                for t in [3000., 2000.]:
                    for p in [1., .1, .01]:
                        inputs.append(("TP", t, p * 101325, [1., of], ()))
        else:
            of = reactants.weight_eq_ratio_to_of_ratio(oxidant, fuel, 1.)
            for rho in [.091864, .0080877, .00066054]:
                inputs.append(("TV", 3000., 1 / rho, [1., of], ()))
    elif example in (3, 4):
        weights, temps = [17., .4, .6], [700., 298.15, 298.15]
        h = reactants.calc_property(cea.ENTHALPY, np.array(weights), np.array(temps))
        if example == 3:
            inputs = [("HP", h, p * 1e5, weights, temps) for p in [100., 10., 1.]]
        else:
            inputs = [("UV", -45.1343 * cea.R, 1 / 14.428, weights, ())]
    elif example == 5:
        weights, temps = [.7206, .1858, .09, .002, .0016], [298.15] * 5
        h = reactants.calc_property(cea.ENTHALPY, np.array(weights), np.array(temps))
        inputs = [("HP", h, p * 1e5, weights, temps) for p in [34.473652, 17.236826, 8.618413, 3.447365, .344737]]
    else:
        weights = reactants.moles_to_weights(np.array([100., 60.]))
        inputs = [("TP", t, .05 * 101325, weights, ()) for t in [1000., 500., 350., 305., 304.3, 304.2, 304., 300., 200., 180., 160.]]
    for i, (mode, first, second, weights, temps) in enumerate(inputs):
        solution = cea.EqSolution(solver)
        solver.solve(solution, getattr(cea, mode), first if mode.startswith("T") else first / cea.R,
                     second / 1e5 if mode.endswith("P") else second, np.array(weights))
        record(f"RP1311-{example}-{i}", names, weights, mode, first, second, solution,
               only=only, omit=omit, trace=trace, temps=temps)

# NASA's generated high-temperature gas fits and both global temperature bounds.
for names, weights in [(["H2", "O2"], [1., 8.]), (["CH4", "O2", "N2"], [1., 4., 13.])]:
    reactants = cea.Mixture(names)
    products = cea.Mixture(names, products_from_reactants=True)
    solver = cea.EqSolver(products, reactants=reactants)
    for t in [6000.1, 10000., 20000., 22000.]:
        solution = cea.EqSolution(solver)
        solver.solve(solution, cea.TP, t, 10., np.array(weights))
        record(f"high-temperature:{names}:{t}", names, weights, "TP", t, 1e6, solution)

args.output.parent.mkdir(parents=True, exist_ok=True)
header = "# Native cea==3.3.4, SI units; generated by scripts/record_equilibrium_parity.py.\n"
header += "# NASA reference 2f79a647e85737742fb08a127da00d1ea1d069da. Relative property tolerance 1e-4, absolute mole-fraction tolerance 3e-5.\n"
header += "# label|reactants|weights|mode|state1|state2|trace|only|omit|reactant temperatures|T;P;rho;H;S;Cp_eq;gamma_s;M|mole fractions\n"
args.output.write_text(header + "\n".join(rows) + "\n")
print(f"Recorded {len(rows)} native states in {args.output}")
