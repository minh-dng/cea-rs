"""Run identical Python calls against NASA 3.3.4 and an installed Rust wheel.

Without NASA, independent stored native values still cover every problem type.
Run: python tests/test_nasa_api.py [--installed /path/to/site-packages]
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("--installed", type=Path)
parser.add_argument("--report", type=Path)
args, remaining = parser.parse_known_args()
if args.installed:
    sys.path.insert(0, str(args.installed.resolve()))

import numpy as np
from cea_rs import nasa as rust
import cea_rs

try:
    import cea as native
except ModuleNotFoundError as exc:
    if exc.name != "cea":
        raise
    native = None

if native is not None and native.__version__ != "3.3.4":
    raise RuntimeError(f"Require cea==3.3.4, got {native.__version__}")

REPORT = {
    "native_version": None if native is None else native.__version__,
    "worst_relative": {},
    "exceptions": [],
}
HERE = Path(__file__).parent
NATIVE_VALUES = json.loads(
    (HERE / "data" / "nasa" / "python-api-3.3.4.json").read_text()
)["cases"]
THERMO = (
    "T",
    "P",
    "density",
    "M",
    "MW",
    "enthalpy",
    "energy",
    "entropy",
    "gibbs_energy",
    "gamma_s",
    "cp_fr",
    "cp_eq",
    "cv_fr",
    "cv_eq",
)
TRANSPORT = ("viscosity", "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq")
PERFORMANCE = (
    "Mach",
    "sonic_velocity",
    "ae_at",
    "c_star",
    "coefficient_of_thrust",
    "Isp",
    "Isp_vacuum",
)
RP1_PRODUCTS = (
    "Ar C CO CO2 H H2 H2O HNO HO2 HNO2 HNO3 N NH NO N2 N2O3 O O2 OH O3".split()
)
RP12_PRODUCTS = "CO CO2 H HNO HO2 H2 H2O H2O2 N NO NO2 N2 N2O O OH O2 HCO NH CH4 NH2 NH3 H2O(L) C(gr)".split()


def mixtures(api, names, only=None, omit=None, ions=False):
    reac = api.Mixture(names, ions=ions)
    prod = (
        api.Mixture(only, ions=ions)
        if only
        else api.Mixture(
            names, products_from_reactants=True, omit=omit or [], ions=ions
        )
    )
    return reac, prod


def equilibrium(api, mode="TP", first=2500.0, second=10.0, *, transport=False):
    reac, prod = mixtures(api, ["H2", "O2"])
    solver = api.EqSolver(prod, reactants=reac, transport=transport)
    soln = api.EqSolution(solver)
    result = solver.solve(soln, getattr(api, mode), first, second, np.array([1.0, 8.0]))
    assert result is None
    return soln, prod, solver


def rp_rocket(api, example=8, n_frz=None, *, transport=False):
    if example in (8, 9, 10):
        names, w, temps = ["H2(L)", "O2(L)"], [1.0, 5.55157], [20.27, 90.17]
        pc, only, ions, trace, insert = 53.3172, None, False, -1.0, []
        options = {"pi_p": [10.0, 100.0, 1000.0], "supar": [25.0, 50.0, 75.0]}
        if example == 8:
            options["subar"] = [1.58]
        elif example == 9:
            options.update(iac=False, ac_at=1.58)
        else:
            options.update(iac=False, mdot=1333.9)
    elif example == 11:
        names, w, temps = ["Li(cr)", "F2(L)"], None, [298.15, 85.02]
        pc, only, ions, trace, insert = (
            api.units.psi_to_bar(1000.0),
            None,
            True,
            -1.0,
            [],
        )
        options = {"pi_p": [68.0457], "subar": [10.0], "supar": [10.0, 20.0, 100.0]}
    elif example == 12:
        names, w, temps = ["CH6N2(L)", "N2O4(L)"], [1.0, 2.5], [298.15, 298.15]
        pc, only, ions, trace, insert = (
            api.units.psi_to_bar(1000.0),
            RP12_PRODUCTS,
            False,
            -1.0,
            [],
        )
        options = {
            "pi_p": [68.0457],
            "supar": [5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0],
        }
    else:
        names, w, temps = (
            ["N2H4(L)", "Be(a)", "H2O2(L)"],
            [0.536, 0.134, 0.33],
            [298.15] * 3,
        )
        pc, only, ions, trace, insert = (
            api.units.psi_to_bar(3000.0),
            None,
            False,
            1e-10,
            ["BeO(L)"],
        )
        options = {"pi_p": [3.0, 10.0, 30.0, 300.0]}
    reac, prod = mixtures(api, names, only=only, ions=ions)
    weights = (
        reac.moles_to_weights(np.array([1.0, 0.5556])) if w is None else np.array(w)
    )
    hc = reac.calc_property(api.ENTHALPY, weights, np.array(temps)) / api.R
    solver = api.RocketSolver(
        prod, reactants=reac, transport=transport, ions=ions, trace=trace, insert=insert
    )
    soln = api.RocketSolution(solver)
    assert solver.solve(soln, weights, pc, hc=hc, n_frz=n_frz, **options) is None
    return soln, prod, solver


def shock(
    api,
    *,
    reflected=False,
    incident_frozen=False,
    reflected_frozen=False,
    speed=1400.0,
    transport=False,
):
    reac, prod = mixtures(api, ["H2", "O2", "Ar"])
    weights = reac.moles_to_weights(np.array([0.05, 0.05, 0.90]))
    solver = api.ShockSolver(prod, reactants=reac, transport=transport)
    soln = api.ShockSolution(solver, reflected=reflected)
    assert (
        solver.solve(
            soln,
            weights,
            300.0,
            0.01333,
            u1=speed,
            reflected=reflected,
            incident_frozen=incident_frozen,
            reflected_frozen=reflected_frozen,
        )
        is None
    )
    return soln, prod, solver


def detonation(api, temperature=298.15, pressure=1.0, *, transport=False):
    reac, prod = mixtures(api, ["H2", "O2"])
    of = reac.chem_eq_ratio_to_of_ratio(np.array([0.0, 1.0]), np.array([1.0, 0.0]), 1.0)
    weights = reac.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), of)
    solver = api.DetonationSolver(prod, reactants=reac, transport=transport)
    soln = api.DetonationSolution(solver)
    assert solver.solve(soln, weights, temperature, pressure) is None
    return soln, prod, solver


def condensed(api, name, mode, first, weights):
    reac, prod = mixtures(api, ["Ar", name])
    solver = api.EqSolver(prod, reactants=reac)
    soln = api.EqSolution(solver)
    solver.solve(soln, getattr(api, mode), first, 1.0, np.array(weights))
    return soln, prod, solver


def air_jet_a(api):
    reac, prod = mixtures(api, ["Air", "Jet-A(L)"])
    weights = np.array([25.0, 1.0])
    hc = reac.calc_property(api.ENTHALPY, weights, np.array([298.15, 298.15])) / api.R
    solver = api.RocketSolver(prod, reactants=reac)
    soln = api.RocketSolution(solver)
    solver.solve(soln, weights, 200.0, hc=hc, supar=80.0)
    return soln, prod, solver


class NasaApiTests(unittest.TestCase):
    def test_derivative_closure_check(self):
        solution, _, solver = equilibrium(rust)
        derivatives = rust.EqDerivatives(solver, solution)
        derivatives.compute_derivatives()
        expected = derivatives.dnj_dw0
        derivatives.compute_derivatives(check_closure_defect=True)
        np.testing.assert_array_equal(derivatives.dnj_dw0, expected)
        core = solver._inner.derivatives("TP", 2500.0, 1e6, [1.0, 8.0])
        rows = {
            name: getattr(core, name)
            for name in (
                "temperature", "total_moles", "enthalpy", "energy",
                "gibbs_energy", "entropy", "species_amounts",
            )
        }
        for defect in (None, float("nan"), 1e-4):
            with self.subTest(defect=defect):
                with patch.object(
                    solver._inner,
                    "derivatives",
                    return_value=SimpleNamespace(**rows, closure_defect=defect),
                ):
                    with self.assertRaisesRegex(cea_rs.ConvergenceError, "closure defect"):
                        derivatives.compute_derivatives(check_closure_defect=True)

    def close(self, group, label, actual, expected, tolerance=1e-4, atol=0):
        a, b = np.asarray(actual), np.asarray(expected)
        self.assertEqual(a.shape, b.shape, label)
        self.assertTrue(np.all(np.isfinite(a)), label)
        relative = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
        meaningful = np.abs(b) > atol
        worst = float(np.max(np.where(meaningful, relative, 0))) if a.size else 0.0
        if worst > REPORT["worst_relative"].get(group, {}).get("difference", -1):
            REPORT["worst_relative"][group] = {"difference": worst, "member": label}
        np.testing.assert_allclose(a, b, rtol=tolerance, atol=atol, err_msg=label)

    def stored(self, key, solution):
        self.assertTrue(solution.converged)
        self.assertEqual(solution.last_error, 0)
        for name, expected in NATIVE_VALUES[key].items():
            self.close(
                "stored NASA " + key,
                name,
                getattr(solution, name),
                expected,
                atol=1e-10 if name in ("enthalpy", "energy") else 0,
            )

    def compare(
        self,
        group,
        actual,
        expected,
        fields=THERMO,
        *,
        frozen=False,
        example12=False,
        reflected_frozen=False,
    ):
        a, ap, solver = actual
        b, bp, _ = expected
        self.assertTrue(a.converged)
        self.assertTrue(b.converged)
        self.assertEqual(a.last_error, 0)
        self.assertEqual(ap.species_names, bp.species_names)
        for name in fields:
            av, bv = getattr(a, name), getattr(b, name)
            self.assertIs(type(av), type(bv), f"{group}.{name} return type")
            if isinstance(av, np.ndarray):
                self.assertEqual(av.dtype, np.float64)
            if frozen and name in ("cv_fr", "cv_eq"):
                # NASA retains a previous Cv or zero after freezing. Check the
                # returned Rust Cv independently; bound the native legacy values.
                self.close(
                    group + " frozen Cv identity",
                    name,
                    a.cv_fr,
                    a.cp_fr - rust.R * a.n / 1000,
                    atol=1e-6,
                )
                eq = np.array([s.equilibrium for s in a._states])
                self.close(group, name + " equilibrium stations", av[eq], bv[eq])
                for i in np.flatnonzero(~eq):
                    retained = np.r_[0.0, bv[: i + 1][eq[: i + 1]], bv[0]]
                    self.assertTrue(
                        np.any(np.isclose(bv[i], retained, rtol=1e-4, atol=1e-12)), name
                    )
                REPORT["exceptions"].append(
                    f"{group}.{name}: native frozen Cv retains a previous value or zero; Rust evaluates current Cv"
                )
                continue
            if (
                reflected_frozen
                and name in ("cp_fr", "cp_eq", "cv_fr", "cv_eq")
                and not (name == "cp_fr" and b.viscosity[-1] > 0)
            ):
                self.close(group, name + " stations 1/2", av[:2], bv[:2], atol=1e-12)
                masses = np.array([a.mass_fractions[n][-1] for n in ap.species_names])
                cp = ap.calc_property(rust.FROZEN_CP, masses, float(a.T[-1])) / 1000
                self.close(group + " station5 Cp", "Cp identity", a.cp_fr[-1], cp)
                self.close(
                    group + " station5 Cv",
                    "Cv identity",
                    a.cv_fr[-1],
                    cp - rust.R * a.n[-1] / 1000,
                )
                REPORT["exceptions"].append(
                    f"{group}.{name}[2]: native reflected-frozen heat capacity is uninitialized; independent fixed-mixture check"
                )
                continue
            self.close(
                group, name, av, bv, atol=1e-10 if name in ("enthalpy", "energy") else 0
            )
        for fractions in ("mole_fractions", "mass_fractions"):
            af, bf = getattr(a, fractions), getattr(b, fractions)
            self.assertEqual(list(af), list(bf))
            for species in af:
                av, bv = np.atleast_1d(af[species]), np.atleast_1d(bf[species])
                mask = (
                    np.maximum(
                        np.atleast_1d(a.mole_fractions[species]),
                        np.atleast_1d(b.mole_fractions[species]),
                    )
                    >= 1e-5
                )
                if example12 and species in ("CO", "OH", "H2"):
                    at75 = np.isclose(a.ae_at, 75.0, rtol=0, atol=1e-8)
                    normal = mask & ~at75
                    self.close(group, f"{fractions}.{species}", av[normal], bv[normal])
                    self.close(
                        "RP12 bounded trace exception",
                        f"{fractions}.{species}[area75]",
                        av[mask & at75],
                        bv[mask & at75],
                        tolerance=1.8e-4,
                    )
                    REPORT["exceptions"].append(
                        f"RP12 {fractions}.{species} at area75: native stops at area74.99733854; bound 1.8e-4"
                    )
                else:
                    self.close(group, f"{fractions}.{species}", av[mask], bv[mask])
        self.assertEqual(
            a.nj.shape,
            (
                (a.num_pts, solver.num_products)
                if isinstance(a, rust.RocketSolution)
                else (len(a.T), solver.num_products)
                if isinstance(a, rust.ShockSolution)
                else (solver.num_products,)
            ),
        )
        gas = ap._gas
        if isinstance(a, rust.ShockSolution):
            with self.assertRaisesRegex(AttributeError, "_get_size"):
                _ = b.nj
            amounts = np.array(
                [
                    b.mass_fractions[n] / r.molecular_weight
                    for n, r in zip(ap.species_names, ap._reactants)
                ]
            ).T
            self.close(
                group,
                "n reconstructed from native fractions",
                a.n,
                amounts[..., gas].sum(axis=-1),
            )
            REPORT["exceptions"].append(
                "NASA ShockSolution.n/nj/ln_nj call missing _get_size; Rust returns n/nj arrays; independent native-fraction check"
            )
        else:
            self.close(group, "n", a.n, b.n)
            amounts = b.nj
        mask = np.array(
            [
                np.maximum(
                    np.atleast_1d(a.mole_fractions[n]),
                    np.atleast_1d(b.mole_fractions[n]),
                )
                >= 1e-5
                for n in ap.species_names
            ]
        ).T
        aa, bb = np.atleast_2d(a.nj), np.atleast_2d(amounts)
        if example12:
            for j, n in enumerate(ap.species_names):
                if n in ("CO", "OH", "H2"):
                    bounded = mask[:, j] & np.isclose(a.ae_at, 75.0, atol=1e-8, rtol=0)
                    self.close(
                        "RP12 bounded trace exception",
                        "nj." + n,
                        aa[bounded, j],
                        bb[bounded, j],
                        tolerance=1.8e-4,
                    )
                    mask[bounded, j] = False
        self.close(group, "nj", aa[mask], bb[mask])
        if all(s.equilibrium for s in a._states):
            self.assertEqual(
                a.ln_nj.shape,
                ((a.num_pts, solver.num_gas) if a._array else (solver.num_gas,)),
            )
            self.close(
                group,
                "ln_nj",
                np.atleast_2d(a.ln_nj)[mask[:, gas]],
                np.atleast_2d(b.ln_nj)[mask[:, gas]],
            )
        self.close(
            group,
            "fraction sum",
            sum(a.mass_fractions.values()),
            np.ones_like(a.T),
            atol=1e-12,
        )

    def test_equilibrium_and_derivatives(self):
        base, _, _ = equilibrium(native or rust)
        cases = {
            "TP": (2500.0, 10.0),
            "HP": (base.enthalpy * 1000 / rust.R, 10.0),
            "SP": (base.entropy * 1000 / rust.R, 10.0),
            "TV": (2500.0, base.volume),
            "UV": (base.energy * 1000 / rust.R, base.volume),
            "SV": (base.entropy * 1000 / rust.R, base.volume),
        }
        for mode, (first, second) in cases.items():
            with self.subTest(mode=mode):
                actual = equilibrium(rust, mode, first, second)
                if native:
                    self.compare(mode, actual, equilibrium(native, mode, first, second))
                else:
                    self.close(mode, "T stored native", actual[0].T, 2500.0)
                    self.close(
                        mode, "H stored native", actual[0].enthalpy, -7414.587537578467
                    )
                d = rust.EqDerivatives(actual[2], actual[0])
                d.compute_derivatives(check_closure_defect=True)
                d.compute_fd(1e-5, central=True)
                if native:
                    oracle = equilibrium(native, mode, first, second)
                    nd = native.EqDerivatives(oracle[2], oracle[0])
                    nd.compute_derivatives()
                    nd.compute_fd(1e-5, central=True)
                    for symbol in ("T", "n", "H", "U", "G", "S"):
                        for variable in ("state1", "state2", "w0"):
                            name = f"d{symbol}_d{variable}"
                            self.close(
                                "derivatives",
                                mode + "." + name,
                                getattr(d, name),
                                getattr(nd, name),
                                atol=1e-7,
                            )
                            fixed = {"HP": "H", "SP": "S", "UV": "U", "SV": "S"}.get(
                                mode
                            )
                            if symbol == fixed and variable != "state1":
                                # Dividing the forward solver's residual by a
                                # 1e-5 step amplifies noise in exact-zero responses.
                                bound = 1e-5 if symbol == "S" else 1e-3
                                for label, obj in (("Rust", d), ("NASA", nd)):
                                    value = getattr(obj, name + "_fd")
                                    self.close(
                                        "FD fixed-constraint noise",
                                        mode + "." + name + " " + label,
                                        value,
                                        np.zeros_like(value),
                                        atol=bound,
                                    )
                                REPORT["exceptions"].append(
                                    f"{mode}.{name}_fd at h=1e-5: exact-zero constraint derivative; both backends bounded absolutely at {bound} in NASA units"
                                )
                            else:
                                self.close(
                                    "native derivative FD",
                                    mode + "." + name + "_fd",
                                    getattr(d, name + "_fd"),
                                    getattr(nd, name + "_fd"),
                                    atol=1e-7,
                                )
                    for variable in ("state1", "state2", "w0"):
                        name = "dnj_d" + variable
                        self.close(
                            "derivatives",
                            mode + "." + name,
                            getattr(d, name),
                            getattr(nd, name),
                            atol=1e-10,
                        )
                self.close(
                    "derivative FD", mode + ".dT_dw0", d.dT_dw0_fd, d.dT_dw0, atol=1e-6
                )

    def test_mixture_helpers_and_units(self):
        for api in (rust, native) if native else (rust,):
            reac = api.Mixture(["CH3OH", "O2"])
            ox, fu = np.array([0.0, 1.0]), np.array([1.0, 0.0])
            r = reac.chem_eq_ratio_to_of_ratio(ox, fu, 0.75)
            phi = reac.weight_eq_ratio_to_of_ratio(ox, fu, 0.75)
            self.close("mixture", "chemical r", r, 2.163754954716528, tolerance=1e-12)
            self.close("mixture", "mass phi", phi, 1.9973122658921798, tolerance=1e-12)
            self.close(
                "mixture",
                "mass/mole round trip",
                reac.moles_to_weights(reac.weights_to_moles(np.array([3.0, 7.0]))),
                [3.0, 7.0],
                tolerance=1e-12,
            )
        if native:
            a, b = rust.Mixture(["H2", "O2"]), native.Mixture(["H2", "O2"])
            for prop in (
                rust.ENTHALPY,
                rust.ENERGY,
                rust.FROZEN_CP,
                rust.FROZEN_CV,
                rust.ENTROPY,
                rust.GIBBS_ENERGY,
                rust.VOLUME,
                rust.DENSITY,
            ):
                for temp in (1000.0, np.array([1000.0, 1200.0])):
                    self.close(
                        "mixture",
                        f"calc_property {prop} {type(temp).__name__}",
                        a.calc_property(
                            prop, np.array([1.0, 8.0]), temp, pressure=10.0
                        ),
                        b.calc_property(
                            prop, np.array([1.0, 8.0]), temp, pressure=10.0
                        ),
                    )
            for name in native.units.__all__:
                self.close(
                    "units",
                    name,
                    getattr(rust.units, name)(np.array([1.25, 10.5])),
                    getattr(native.units, name)(np.array([1.25, 10.5])),
                    tolerance=1e-14,
                )
            for name in native.__all__:
                if isinstance(getattr(native, name), int):
                    self.assertEqual(getattr(rust, name), getattr(native, name), name)
        for api in (rust, native) if native else (rust,):
            custom = api.Reactant(
                "binder",
                formula={"C": 1.0, "H": 2.0},
                enthalpy=-100.0,
                enthalpy_units="kJ/mol",
                temperature=298.15,
            )
            mix = api.Mixture([custom, "O2"])
            self.close(
                "custom",
                "molar enthalpy",
                mix.calc_property(api.ENTHALPY, np.array([1.0, 0.0]), 298.15),
                -100000000 / 14.02658,
            )

    def test_rp_equilibrium_and_condensed_exceptions(self):
        def rp(api, mode, temperature, second, names, weights, only=None):
            reac, prod = mixtures(api, names, only=only)
            solver = api.EqSolver(prod, reactants=reac, transport=mode == "TV")
            soln = api.EqSolution(solver)
            solver.solve(
                soln, getattr(api, mode), temperature, second, np.array(weights)
            )
            return soln, prod, solver

        # RP1/2 explicit gas pool and RP14's water condensation schedule.
        for mode, t, second in (
            ("TP", 3000.0, rust.units.atm_to_bar(1.0)),
            ("TV", 3000.0, 1 / 0.091864),
        ):
            reac = rust.Mixture(["H2", "Air"])
            weights = reac.of_ratio_to_weights(
                np.array([0.0, 1.0]),
                np.array([1.0, 0.0]),
                reac.weight_eq_ratio_to_of_ratio(
                    np.array([0.0, 1.0]), np.array([1.0, 0.0]), 1.0
                ),
            )
            a = rp(rust, mode, t, second, ["H2", "Air"], weights, RP1_PRODUCTS)
            self.stored(f"rp-{mode}", a[0])
            if native:
                self.compare(
                    mode,
                    a,
                    rp(native, mode, t, second, ["H2", "Air"], weights, RP1_PRODUCTS),
                    THERMO + (TRANSPORT if mode == "TV" else ()),
                )
        for t in (1000.0, 500.0, 350.0, 305.0, 304.3, 304.2, 304.0, 300.0):
            weights = rust.Mixture(["H2(L)", "O2(L)"]).moles_to_weights(
                np.array([100.0, 60.0])
            )
            a = rp(
                rust, "TP", t, rust.units.atm_to_bar(0.05), ["H2(L)", "O2(L)"], weights
            )
            self.stored(f"rp14-{t}", a[0])
            if native:
                self.compare(
                    "condensed TP",
                    a,
                    rp(
                        native,
                        "TP",
                        t,
                        rust.units.atm_to_bar(0.05),
                        ["H2(L)", "O2(L)"],
                        weights,
                    ),
                )
        for name, mode, first, w, temperature, phase, fraction in (
            (
                "Hg",
                "HP",
                -15551.803948663037 / rust.R,
                [1.0, 20.0],
                234.28999933688,
                "Hg(L)",
                0.479137243881,
            ),
            ("H2O", "TP", 273.15, [1.0, 1.0], 273.15, "H2O(L)", 0.498612290857),
            ("NaCN", "TP", 288.5, [1.0, 1.0], 288.5, "NaCN(III)", 0.499999999526),
        ):
            a = condensed(rust, name, mode, first, w)[0]
            self.assertTrue(a.converged)
            self.close("condensed stored", name + " T", a.T, temperature)
            self.close(
                "condensed stored", name + " phase", a.mass_fractions[phase], fraction
            )
            if native:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    b = condensed(native, name, mode, first, w)[0]
                self.assertFalse(b.converged)
                self.assertEqual(b.last_error, native.NOT_CONVERGED)
            REPORT["exceptions"].append(
                f"condensed {name} {mode} {first}: native error8; Rust matches stored independent T/phase checks"
            )
        a = condensed(rust, "Br2", "TP", 262.0, [1.0, 20.0])[0]
        self.close(
            "condensed stored",
            "Br2 solid mass",
            a.mass_fractions["Br2(cr)"],
            0.943263839016,
        )
        if native:
            b = condensed(native, "Br2", "TP", 262.0, [1.0, 20.0])[0]
            self.assertTrue(b.converged)
            self.assertEqual(b.mass_fractions["Br2(cr)"], 0)
            self.assertGreater(b.mole_fractions["Br2"], 10 * a.mole_fractions["Br2"])
        REPORT["exceptions"].append(
            "Br2 TP262K: native reversed interval omits solid; Rust solid mass0.943263839016 and lower saturated vapour"
        )

    def test_air_jet_a_false_positive(self):
        a, _, _ = air_jet_a(rust)
        self.assertTrue(a.converged)
        self.close("Air Jet-A stored", "exit T", a.T[-1], 273.121944437)
        self.close("Air Jet-A stored", "exit P", a.P[-1], 0.09971217843957)
        self.close("Air Jet-A invariants", "entropy", a.entropy[-1], a.entropy[0])
        self.close(
            "Air Jet-A invariants",
            "energy",
            a.enthalpy[-1] + a.Isp[-1] ** 2 / 2000,
            a.enthalpy[0],
            atol=1e-5,
        )
        self.close(
            "Air Jet-A invariants",
            "mass flux",
            a.density[-1] * a.Isp[-1] * 80,
            a.density[1] * a.Isp[1],
        )
        if native:
            b, _, _ = air_jet_a(native)
            self.assertTrue(b.converged)
            residual = b.enthalpy[-1] + b.Isp[-1] ** 2 / 2000 - b.enthalpy[0]
            self.assertGreater(residual, 14.0)
            self.assertLess(residual, 15.0)
            flux = b.density[-1] * b.Isp[-1] * 80 / (b.density[1] * b.Isp[1]) - 1
            self.assertGreater(flux, 0.18)
            self.assertLess(flux, 0.19)
        REPORT["exceptions"].append(
            "Air/Jet-A 20MPa area80: NASA success has 14–15kJ/kg energy error and 18–19% mass-flux error; Rust independently closes both"
        )

    def test_rocket_schedules(self):
        for example, nf, transport in (
            (8, None, False),
            (8, 1, True),
            (8, 2, True),
            (9, None, False),
            (9, 4, False),
            (9, 5, False),
            (10, None, False),
            (11, None, True),
            (12, None, False),
            (12, 2, False),
            (13, None, False),
        ):
            with self.subTest(example=example, n_frz=nf):
                actual = rp_rocket(rust, example, nf, transport=transport)
                self.stored(f"rocket-{example}-{nf or 0}", actual[0])
                if native:
                    self.compare(
                        ("FAC frozen" if nf else "FAC equilibrium")
                        if example in (9, 10)
                        else "IAC frozen"
                        if nf
                        else "IAC equilibrium",
                        actual,
                        rp_rocket(native, example, nf, transport=transport),
                        THERMO + PERFORMANCE + (TRANSPORT if transport else ()),
                        frozen=nf is not None,
                        example12=example == 12,
                    )
        # Native's early FAC freeze fails and retains seven stations; accept
        # only this exact failure and verify the complete Rust solution's energy.
        a = rp_rocket(rust, 9, 3)[0]
        self.assertTrue(a.converged)
        self.assertEqual(a.num_pts, 10)
        self.close(
            "FAC early freeze invariant",
            "stagnation enthalpy",
            a.enthalpy[2:] + a.Isp[2:] ** 2 / 2000,
            np.full(a.num_pts - 2, a.enthalpy[0]),
            atol=1e-6,
        )
        if native:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                b = rp_rocket(native, 9, 3)[0]
            self.assertFalse(b.converged)
            self.assertEqual(b.last_error, native.NOT_CONVERGED)
            self.assertEqual(b.num_pts, 7)
        REPORT["exceptions"].append(
            "RP9 n_frz=3: native error8 stops after seven stations; Rust returns ten verified stations with energy closure"
        )

    def test_shock_and_detonation(self):
        for reflected, incident, reflection, transport in product(
            (False, True), repeat=4
        ):
            with self.subTest(
                reflected=reflected,
                incident=incident,
                reflection=reflection,
                transport=transport,
            ):
                kw = dict(
                    reflected=reflected,
                    incident_frozen=incident,
                    reflected_frozen=reflection,
                    transport=transport,
                )
                a = shock(rust, **kw)
                self.stored(
                    f"shock-{int(reflected)}{int(incident)}{int(reflection)}{int(transport)}",
                    a[0],
                )
                if native:
                    self.compare(
                        "shock",
                        a,
                        shock(native, **kw),
                        THERMO
                        + TRANSPORT
                        + (
                            "velocity",
                            "Mach",
                            "sonic_velocity",
                            "rho12",
                            "rho52",
                            "P21",
                            "P52",
                            "T21",
                            "T52",
                            "M21",
                            "M52",
                            "v2",
                            "u5_p_v2",
                        ),
                        reflected_frozen=reflected and reflection,
                    )
        for t, p in ((298.15, 1.0), (298.15, 20.0), (500.0, 1.0), (500.0, 20.0)):
            a = detonation(rust, t, p, transport=True)
            self.stored(f"detonation-{t}-{p}", a[0])
            if native:
                self.compare(
                    "detonation",
                    a,
                    detonation(native, t, p, transport=True),
                    THERMO
                    + TRANSPORT
                    + (
                        "P1",
                        "T1",
                        "H1",
                        "M1",
                        "gamma1",
                        "sonic_velocity1",
                        "Mach",
                        "velocity",
                        "sonic_velocity",
                        "P_P1",
                        "T_T1",
                        "M_M1",
                        "rho_rho1",
                    ),
                )

    def test_failures_reuse_and_unsupported(self):
        s, _, solver = equilibrium(rust)
        for w in ([], [0, 0], [-1, 2], [float("nan"), 1]):
            with self.assertRaises(ValueError):
                solver.solve(s, rust.TP, 3000.0, 10.0, np.array(w))
            self.assertFalse(s.converged)
            with self.assertRaisesRegex(RuntimeError, "no verified solution"):
                _ = s.T
        reac, prod = mixtures(rust, ["H2", "O2"])
        for kwargs in ({"smooth_truncation": True}, {"truncation_width": 0.5}):
            with self.assertRaisesRegex(NotImplementedError, "truncation"):
                rust.EqSolver(prod, reactants=reac, **kwargs)
        for call, feature in (
            (lambda: rust.init(path="anything"), "database"),
            (lambda: rust.init_thermo("thermo.lib"), "thermodynamic"),
            (lambda: rust.init_trans("trans.lib"), "transport"),
            (lambda: rust.EqSolution(solver, T_init=2500), "T_init"),
        ):
            with self.assertRaisesRegex(NotImplementedError, feature):
                call()
        ds = rust.DetonationSolver(prod, reactants=reac)
        with self.assertRaisesRegex(NotImplementedError, "Frozen detonation"):
            ds.solve(
                rust.DetonationSolution(ds),
                np.array([1.0, 8.0]),
                300.0,
                1.0,
                frozen=True,
            )
        with self.assertRaises(cea_rs.ConvergenceError):
            shock(rust, speed=1000.0)
        if native:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                n = shock(native, speed=1000.0)[0]
            self.assertFalse(n.converged)
            self.assertEqual(n.last_error, native.LAST_VALID_SOLUTION)
        REPORT["exceptions"].append(
            "RP7 1000 m/s equilibrium: native LAST_VALID_SOLUTION=9; Rust raises and exposes no state"
        )

    def test_helpers_defaults_and_namespace(self):
        from cea_rs.nasa.lib import libcea
        from cea_rs.nasa import matlab

        self.assertIs(libcea.EqSolver, rust.EqSolver)
        self.assertEqual(rust.lib_version(), cea_rs.__version__)
        self.assertEqual(
            (
                rust.lib_version_major(),
                rust.lib_version_minor(),
                rust.lib_version_patch(),
            ),
            (0, 1, 0),
        )
        self.assertIn("nasa_python_api", cea_rs.capabilities())
        self.assertTrue(rust.is_initialized())
        rust.init(trans=False)
        rust.set_log_level(rust.LOG_NONE)
        with self.assertRaisesRegex(NotImplementedError, "logging"):
            rust.set_log_level(rust.LOG_DEBUG)
        with self.assertRaisesRegex(NotImplementedError, "MATLAB"):
            _ = matlab.rocket_solve
        with self.assertRaisesRegex(NotImplementedError, "fitted named reactant"):
            rust.Mixture([rust.Reactant("H2", enthalpy=1.0, enthalpy_units="kJ/mol")])

        def helper(api, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                return api.eq_solve(
                    api.TP,
                    ["H2", "O2"],
                    T=2500.123456789,
                    P=10.0,
                    fuel_amounts=[1.0, 0.0],
                    oxid_amounts=[0.0, 1.0],
                    only=["H", "H2", "H2O", "OH", "O", "O2"],
                    **kwargs,
                )

        for kwargs in (
            {"of_ratio": 8.0},
            {"phi": 1.0},
            {"r_eq": 0.75},
            {"pct_fuel": 20.0},
            {"moles": True},
        ):
            a = helper(rust, **kwargs)
            self.assertTrue(a.converged)
            self.assertEqual(a.T, 2500.123456789)
            if native:
                b = helper(native, **kwargs)
                self.close("eq_solve", "enthalpy", a.enthalpy, b.enthalpy)

        def defaults(api):
            products = api.Mixture(["H", "H2", "H2O", "OH", "O", "O2"])
            solver = api.EqSolver(products)
            x = api.EqSolution(solver)
            solver.solve(
                x, api.TP, 2500.0, 10.0, np.array([0.0, 1.0, 0.0, 0.0, 0.0, 8.0])
            )
            return x, products, solver

        a = defaults(rust)
        if native:
            self.compare("default reactants", a, defaults(native))
        old = a[0].nj.copy()
        a[0].nj[:] = 0
        np.testing.assert_array_equal(a[0].nj, old)
        with self.assertRaisesRegex(NotImplementedError, "frozen stations"):
            _ = shock(rust)[0].ln_nj

        # omit acts only on generated products, including absent-name no-ops.
        generated = rust.Mixture(["H2", "O2"], True, omit=["H2O(L)", "not-in-database"])
        self.assertNotIn("H2O(L)", generated.species_names)
        explicit = rust.Mixture(["H2", "O2"], omit=["H2"])
        self.assertEqual(explicit.species_names, ["H2", "O2"])
        if native:
            self.assertEqual(
                generated.species_names,
                native.Mixture(
                    ["H2", "O2"], True, omit=["H2O(L)", "not-in-database"]
                ).species_names,
            )


if __name__ == "__main__":
    result = unittest.main(
        argv=[sys.argv[0], *remaining], verbosity=2, exit=False
    ).result
    REPORT["exceptions"] = sorted(set(REPORT["exceptions"]))
    REPORT["tests"] = {
        "run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "skipped": len(result.skipped),
    }
    if args.report:
        args.report.write_text(json.dumps(REPORT, indent=2) + "\n")
    print(json.dumps(REPORT, indent=2))
    sys.exit(not result.wasSuccessful())
