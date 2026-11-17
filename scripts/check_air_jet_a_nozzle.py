"""Audit the Air/Jet-A water-phase exit using independent NASA polynomial sums.

Run with separately installed cea==3.3.4, numpy, and the release cea_rs wheel.
The source checkout supplies only thermo.inp; no native result seeds Rust.
The JSON also records 90 neighboring Rust cases for before/after comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import warnings
from functools import cache
from pathlib import Path

import cea
import cea_rs
import numpy as np

R = 8314.51
ROOT = Path(__file__).resolve().parents[1]
THERMO = ROOT / "vendor/thermo-species/data/thermo.inp"
LINES = THERMO.read_text().splitlines()
WEIGHTS = np.array([25.0, 1.0])


@cache
def species(name: str) -> dict:
    """Read the fixed-width NASA record used at these temperatures (<6000 K)."""
    i = next(i for i, line in enumerate(LINES) if line.split()[:1] == [name])
    meta = LINES[i + 1]
    formula = {
        meta[k:k + 2].strip(): float(meta[k + 2:k + 8])
        for k in range(10, 50, 8) if meta[k:k + 2].strip()
    }
    if name == "Air":
        # NASA's database compiler restores the rounded input composition.
        formula = {"N": 1.56168, "O": 0.419590, "AR": 0.009365, "C": 0.000319}
    intervals = []
    for j in range(int(meta[:2])):
        header, a, b = LINES[i + 2 + 3 * j:i + 5 + 3 * j]
        coefficients = [float(a[k:k + 16].replace("D", "e")) for k in range(0, 80, 16)]
        coefficients += [float(b[k:k + 16].replace("D", "e")) for k in [0, 16, 48, 64]]
        intervals.append((float(header[:11]), coefficients))
    return {"mw": float(meta[52:65]), "phase": int(meta[50:52]),
            "formula": formula, "intervals": intervals}


def thermo(name: str, t: float) -> tuple[float, float]:
    """Return H [J/kmol] and standard S [J/(kmol K)], without solver calls."""
    assert 160 <= t < 6000
    intervals = species(name)["intervals"]
    a = intervals[0][1]
    for lower, coefficients in intervals:
        if t > lower:
            a = coefficients
    h = (-a[0] / t + a[1] * math.log(t) + a[2] * t + a[3] * t**2 / 2
         + a[4] * t**3 / 3 + a[5] * t**4 / 4 + a[6] * t**5 / 5 + a[7])
    s = (-a[0] / (2 * t**2) - a[1] / t + a[2] * math.log(t) + a[3] * t
         + a[4] * t**2 / 2 + a[5] * t**3 / 3 + a[6] * t**4 / 4 + a[8])
    return h * R, s * R


def composition_properties(row: dict) -> dict:
    t, p = row["temperature_k"], row["pressure_pa"]
    amounts = {name: y / species(name)["mw"] for name, y in row["mass_fractions"].items() if y > 0}
    ng = sum(n for name, n in amounts.items() if species(name)["phase"] == 0)
    h, s, elements = 0.0, 0.0, {}
    for name, n in amounts.items():
        data = species(name)
        hj, sj = thermo(name, t)
        mixing = R * math.log(n / ng * p / 1e5) if data["phase"] == 0 else 0.0
        h += n * hj
        s += n * (sj - mixing)
        for element, count in data["formula"].items():
            elements[element] = elements.get(element, 0.0) + n * count
    return {"enthalpy_j_kg": h, "entropy_j_kg_k": s,
            "density_kg_m3": p / (ng * R * t), "elements_kmol_kg": elements,
            "water_partial_pressure_pa": p * amounts.get("H2O", 0.0) / ng}


def checks(rows: list[dict]) -> dict:
    chamber, throat, exit = rows
    properties = composition_properties(exit)
    feed = {}
    for name, w in zip(["Air", "Jet-A(L)"], WEIGHTS / WEIGHTS.sum()):
        for element, count in species(name)["formula"].items():
            feed[element] = feed.get(element, 0.0) + float(w) / species(name)["mw"] * count
    flux = throat["density_kg_m3"] * throat["velocity_m_s"]
    u = exit["velocity_m_s"]
    t = exit["temperature_k"]
    gibbs = {name: thermo(name, t)[0] - t * thermo(name, t)[1]
             for name in ["H2O", "H2O(cr)", "H2O(L)"]}
    saturation = {name: 1e5 * math.exp((gibbs[name] - gibbs["H2O"]) / (R * t))
                  for name in ["H2O(cr)", "H2O(L)"]}
    return properties | {
        "element_relative_errors": {e: properties["elements_kmol_kg"][e] / b - 1 for e, b in feed.items()},
        "entropy_difference_j_kg_k": properties["entropy_j_kg_k"] - chamber["entropy_j_kg_k"],
        "energy_residual_j_kg": properties["enthalpy_j_kg"] + u**2 / 2 - chamber["enthalpy_j_kg"],
        "reported_energy_residual_j_kg": exit["enthalpy_j_kg"] + u**2 / 2 - chamber["enthalpy_j_kg"],
        "throat_mass_flux_kg_m2_s": flux,
        "exit_mass_flux_kg_m2_s": properties["density_kg_m3"] * u,
        "area_from_mass_flux": flux / (properties["density_kg_m3"] * u),
        "area80_relative_residual": properties["density_kg_m3"] * u * 80 / flux - 1,
        "reported_density_area80_relative_residual": exit["density_kg_m3"] * u * 80 / flux - 1,
        "saturation_pressure_pa": saturation,
        "saturation_relative_residual": {name: properties["water_partial_pressure_pa"] / p - 1
                                         for name, p in saturation.items()},
        "g_liquid_minus_ice_j_mol": (gibbs["H2O(L)"] - gibbs["H2O(cr)"]) / 1000,
    }


def native_rows(result: cea.RocketSolution) -> list[dict]:
    fields = {"temperature_k": ("T", 1), "pressure_pa": ("P", 1e5),
              "enthalpy_j_kg": ("enthalpy", 1e3), "entropy_j_kg_k": ("entropy", 1e3),
              "energy_j_kg": ("energy", 1e3), "cp_equilibrium_j_kg_k": ("cp_eq", 1e3),
              "cp_frozen_j_kg_k": ("cp_fr", 1e3),
              "density_kg_m3": ("density", 1), "gamma_s": ("gamma_s", 1),
              "molecular_weight_kg_kmol": ("M", 1), "velocity_m_s": ("Isp", 1),
              "sound_speed_m_s": ("sonic_velocity", 1), "mach": ("Mach", 1),
              "area_ratio": ("ae_at", 1), "vacuum_isp_m_s": ("Isp_vacuum", 1),
              "cstar_m_s": ("c_star", 1), "thrust_coefficient": ("coefficient_of_thrust", 1)}
    return [{key: float(getattr(result, prop)[i]) * scale for key, (prop, scale) in fields.items()}
            | {"mass_fractions": {name: float(y[i]) for name, y in result.mass_fractions.items() if y[i] > 0}}
            | {"mole_fractions": {name: float(x[i]) for name, x in result.mole_fractions.items()}}
            for i in range(result.num_pts)]


def port_rows(result) -> list[dict]:
    fields = ["temperature_k", "pressure_pa", "enthalpy_j_kg", "entropy_j_kg_k",
              "density_kg_m3", "gamma_s", "molecular_weight_kg_kmol", "energy_j_kg",
              "cp_equilibrium_j_kg_k", "cp_frozen_j_kg_k", "mole_fractions"]
    return [{key: getattr(st.state, key) for key in fields}
            | {key: getattr(st, key) for key in ["velocity_m_s", "sound_speed_m_s", "mach", "area_ratio",
                                               "cstar_m_s", "thrust_coefficient"]}
            | {"vacuum_isp_m_s": st.isp_vacuum_m_s,
               "mass_fractions": {name: y for name, y in st.state.mass_fractions.items() if y > 0}}
            for st in result.stations]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
    assert cea.__version__ == "3.3.4"
    reactants = cea.Mixture(["Air", "Jet-A(L)"])
    products = cea.Mixture(["Air", "Jet-A(L)"], products_from_reactants=True)
    rocket_solver = cea.RocketSolver(products, reactants=reactants, transport=False)
    native = cea.RocketSolution(rocket_solver)
    h0 = reactants.calc_property(cea.ENTHALPY, WEIGHTS, np.array([298.15] * 2))
    rocket_solver.solve(native, WEIGHTS, 200.0, supar=[80.0], hc=h0 / R)
    rows = native_rows(native)
    # Calibrate the independent sums against the well-converged hot chamber.
    chamber = composition_properties(rows[0])
    assert abs(chamber["enthalpy_j_kg"] - rows[0]["enthalpy_j_kg"]) < 0.01
    assert abs(chamber["entropy_j_kg_k"] - rows[0]["entropy_j_kg_k"]) < 0.001
    report = {"native_version": cea.__version__, "thermo_sha256": hashlib.sha256(THERMO.read_bytes()).hexdigest(),
              "inputs": {"reactants": ["Air", "Jet-A(L)"], "weights": WEIGHTS.tolist(),
                         "reactant_temperatures_k": [298.15] * 2, "chamber_pressure_pa": 20e6,
                         "supersonic_area_ratio": 80, "condensed": True, "ions": False,
                         "reactant_enthalpy_j_kg": h0},
              "native": {"converged": bool(native.converged), "last_error": native.last_error,
                         "stations": rows, "checks": checks(rows)}}
    eq = cea.EqSolver(products, reactants=reactants)
    report["native_sp"] = []
    for entropy in [rows[-1]["entropy_j_kg_k"], report["native"]["checks"]["entropy_j_kg_k"]]:
        sol = cea.EqSolution(eq)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            eq.solve(sol, cea.SP, entropy / R, rows[-1]["pressure_pa"] / 1e5, WEIGHTS)
        row = {"assigned_entropy_j_kg_k": entropy, "pressure_pa": rows[-1]["pressure_pa"],
               "converged": bool(sol.converged), "last_error": sol.last_error,
               "warnings": [str(w.message) for w in caught]}
        if sol.converged:
            row.update(temperature_k=sol.T, enthalpy_j_kg=sol.enthalpy * 1e3)
        report["native_sp"].append(row)
    solver = cea_rs.Solver(["Air", "Jet-A(L)"])
    try:
        rows = port_rows(solver.rocket(WEIGHTS.tolist(), [298.15] * 2, 20e6, 80))
        report["port"] = {"stations": rows, "checks": checks(rows)}
        report["port_sp"] = []
        for pressure in [report["native"]["stations"][-1]["pressure_pa"], rows[-1]["pressure_pa"]]:
            sol = solver.solve("SP", rows[0]["entropy_j_kg_k"], pressure, WEIGHTS.tolist())
            report["port_sp"].append({"pressure_pa": pressure, "temperature_k": sol.temperature_k,
                                      "entropy_j_kg_k": sol.entropy_j_kg_k, "enthalpy_j_kg": sol.enthalpy_j_kg,
                                      "water_mass_fractions": {name: sol.mass_fractions[name]
                                                               for name in ["H2O", "H2O(cr)", "H2O(L)"]}})
    except cea_rs.ConvergenceError as error:
        report["port"] = {"error": str(error)}
    report["neighbors"] = []
    for of, temp in [(25, 298.15), (15, 298.15), (15, 300.0)]:
        for pressure in [0.7e6, 3e6, 7e6, 10e6, 20e6]:
            for area in [15, 20, 50, 80, 100, 200]:
                row = {"of": of, "reactant_temperature_k": temp, "chamber_pressure_pa": pressure, "area_ratio": area}
                try:
                    result = solver.rocket([of, 1], [temp] * 2, pressure, area)
                    values = port_rows(result)
                    row["stations_sha256"] = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()
                except cea_rs.ConvergenceError as error:
                    row["error"] = str(error)
                report["neighbors"].append(row)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {args.output}; neighboring successes: {sum('error' not in row for row in report['neighbors'])}/90")


if __name__ == "__main__":
    main()
