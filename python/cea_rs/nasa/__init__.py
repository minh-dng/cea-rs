"""NASA cea 3.3.4 calling conventions over the Rust solver.

Use ``from cea_rs import nasa as cea``. This module uses bar, K and NASA's
kJ-based solution properties; the parent cea_rs API continues to use SI.
See docs/nasa-python-api.md for the member inventory and bounded differences.
"""

from __future__ import annotations

import warnings as _warnings

import numpy as _np

import cea_rs as _rs
from cea_rs._cea_rs import _nasa_species

from . import constants, units
from .constants import *

__version__ = _rs.__version__
_DATABASE = {name: (mw, phase, formula) for name, mw, phase, formula in _nasa_species()}
_ATOMIC_MASSES = {
    next(iter(f)): mw
    for mw, phase, f in _DATABASE.values()
    if phase == 0 and len(f) == 1 and next(iter(f.values())) == 1
}
_MODES = ("TP", "HP", "SP", "TV", "UV", "SV")


def lib_version():
    return __version__


def lib_version_major():
    return int(__version__.split(".")[0])


def lib_version_minor():
    return int(__version__.split(".")[1])


def lib_version_patch():
    return int(__version__.split(".")[2].split("a")[0])


def init(path=None, trans=True):
    """Use the embedded databases; custom database paths are unsupported."""
    if path is not None:
        raise NotImplementedError("init(path): custom NASA binary databases")


def init_thermo(thermofile):
    raise NotImplementedError("init_thermo: custom thermodynamic databases")


def init_trans(transfile):
    raise NotImplementedError("init_trans: custom transport databases")


def is_initialized():
    return True


def set_log_level(level):
    if level not in (
        LOG_NONE,
        LOG_DEBUG,
        LOG_INFO,
        LOG_WARNING,
        LOG_ERROR,
        LOG_CRITICAL,
    ):
        raise ValueError("Unknown NASA log level")
    if level in (LOG_DEBUG, LOG_INFO):
        raise NotImplementedError("NASA debug/info logging")
    # Rust has no native logger; Python failures still raise their typed errors.


def _vector(values, length, name="weights", *, positive=False):
    a = _np.asarray(values, dtype=_np.float64)
    if a.ndim != 1 or len(a) != length:
        raise ValueError(
            f"{name} must have {length} entries in a one-dimensional array"
        )
    if not _np.all(_np.isfinite(a)) or _np.any(a < 0) or not _np.isfinite(a.sum()):
        raise ValueError(f"{name} must be finite and nonnegative")
    if positive and a.sum() <= 0:
        raise ValueError(f"{name} must have a positive sum")
    return a


def _scalar32(value):
    # NASA's Cython float arguments round before entering the double-precision core.
    return float(_np.float32(value))


def _state_inputs(eq_type, state1, state2):
    if eq_type not in range(6):
        raise NotImplementedError(f"equilibrium type {eq_type}")
    mode = _MODES[eq_type]
    return (
        mode,
        float(state1) * (1 if mode.startswith("T") else R),
        float(state2) * (1e5 if mode.endswith("P") else 1),
    )


class Reactant:
    """Named or custom reference reactant, with NASA's explicit enthalpy units."""

    def __init__(
        self,
        name,
        formula=None,
        molecular_weight=None,
        enthalpy=None,
        enthalpy_units=None,
        temperature=None,
    ):
        if not isinstance(name, str) or not name.strip() or "\0" in name:
            raise TypeError("Reactant name must be a non-empty str without NUL")
        if formula is not None:
            if not isinstance(formula, dict) or not formula:
                raise ValueError("Reactant formula must be a non-empty dict")
            if any(
                not isinstance(e, str) or not e or not _np.isfinite(float(n))
                for e, n in formula.items()
            ):
                raise ValueError(
                    "Reactant formula must contain element names and finite coefficients"
                )
        for value in (molecular_weight, enthalpy, temperature):
            if value is not None and not _np.isfinite(float(value)):
                raise ValueError("Reactant values must be finite")
        if enthalpy is not None and enthalpy_units is None:
            raise ValueError("Reactant enthalpy requires explicit enthalpy_units")
        if enthalpy_units is not None:
            text = (
                enthalpy_units.decode()
                if isinstance(enthalpy_units, bytes)
                else enthalpy_units
            )
            text = (
                "".join(text.split())
                .lower()
                .replace("/mol", "/mole")
                .replace("/molee", "/mole")
            )
            if text not in {
                f"{e}/{b}" for e in ("j", "kj", "cal", "kcal") for b in ("kg", "mole")
            }:
                raise ValueError(
                    f"Unsupported Reactant enthalpy_units: {enthalpy_units}"
                )
            enthalpy_units = text
            if (
                enthalpy is not None
                and text.endswith("/kg")
                and molecular_weight is None
            ):
                raise ValueError("Weight-based enthalpy requires molecular_weight")
        self.name = name
        self.formula = formula
        self.molecular_weight = molecular_weight
        self.enthalpy = enthalpy
        self.enthalpy_units = enthalpy_units
        self.temperature = temperature

    def get_valid_temperature_range(self):
        """Query the database name, independently of custom overrides (K)."""
        return _rs.Reactant(self.name).get_valid_temperature_range()

    def _native(self):
        if (
            self.formula is None
            and self.enthalpy is None
            and self.molecular_weight is None
        ):
            return _rs.Reactant(self.name)
        try:
            database = _rs.Reactant(self.name)
        except ValueError:
            database = None
        if database is not None:
            lo, hi = database.get_valid_temperature_range()
            if database.thermo((lo + hi) / 2)[0] is not None:
                raise NotImplementedError("overrides of fitted named reactant data")
        named = None
        if self.formula is None:
            named = _rs.Reactant(self.name)
        formula = self.formula if self.formula is not None else named._nasa_formula
        mw = self.molecular_weight
        if mw is None:
            try:
                mw = sum(_ATOMIC_MASSES[e.upper()] * n for e, n in formula.items())
            except KeyError as exc:
                raise ValueError(f"Unknown element {exc.args[0]}") from exc
        h = 0.0 if self.enthalpy is None else float(self.enthalpy)
        if self.enthalpy is not None:
            energy, basis = self.enthalpy_units.split("/")
            h *= {"j": 1, "kj": 1000, "cal": 4.184, "kcal": 4184}[energy]
            if basis == "kg":
                h *= mw / 1000
        return _rs.Reactant(
            self.name,
            formula=formula,
            molecular_weight=float(mw),
            enthalpy_j_mol=h,
            reference_temperature_k=298.15
            if self.temperature is None
            else float(self.temperature),
        )


class Mixture:
    """Species, product selection and fixed-composition mixture helpers."""

    def __init__(
        self, species=None, products_from_reactants=False, omit=[], ions=False
    ):
        if not isinstance(species, list) or not species:
            raise TypeError("Mixture species must be a non-empty list")
        if not all(isinstance(s, (str, Reactant)) for s in species):
            raise TypeError("Mixture species entries must be str or Reactant")
        if omit is not None and not isinstance(omit, list):
            raise TypeError("Mixture omit must be a list")
        self._ions = bool(ions)
        self._reactants = [
            Reactant(s)._native() if isinstance(s, str) else s._native()
            for s in species
        ]
        if products_from_reactants:
            elements = set().union(*(r._nasa_formula for r in self._reactants))
            if ions:
                elements.add("E")
            self._reactants = [
                _rs.Reactant(n)
                for n, (_, _, f) in _DATABASE.items()
                if set(f) <= elements and n not in (omit or [])
            ]
        self._names = [r.name for r in self._reactants]
        self._mw = _np.array([r.molecular_weight for r in self._reactants])
        self._gas = _np.array([r._nasa_phase == 0 for r in self._reactants])
        self.num_species = len(self._names)
        if not self.num_species:
            raise ValueError("Mixture has no species")
        self._helper = _rs.Solver(self._reactants, ions=ions)

    @property
    def species_names(self):
        return self._names.copy()

    def moles_to_weights(self, moles):
        return _np.array(
            self._helper.moles_to_weights(_vector(moles, self.num_species).tolist())
        )

    def weights_to_moles(self, weights):
        return _np.array(
            self._helper.weights_to_moles(_vector(weights, self.num_species).tolist())
        )

    def chem_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, chem_eq_ratio):
        return self._helper.chem_eq_ratio_to_of_ratio(
            _vector(oxidant_weights, self.num_species).tolist(),
            _vector(fuel_weights, self.num_species).tolist(),
            _scalar32(chem_eq_ratio),
        )

    def weight_eq_ratio_to_of_ratio(
        self, oxidant_weights, fuel_weights, weight_eq_ratio
    ):
        return self._helper.weight_eq_ratio_to_of_ratio(
            _vector(oxidant_weights, self.num_species).tolist(),
            _vector(fuel_weights, self.num_species).tolist(),
            _scalar32(weight_eq_ratio),
        )

    def of_ratio_to_weights(self, oxidant_weights, fuel_weights, of_ratio):
        return _np.array(
            self._helper.of_ratio_to_weights(
                _vector(oxidant_weights, self.num_species).tolist(),
                _vector(fuel_weights, self.num_species).tolist(),
                _scalar32(of_ratio),
            )
        )

    def calc_property(self, prop_type, weights, temperature, pressure=None):
        """Fixed composition, retaining 3.3.4's scalar/array unit conventions."""
        w = _vector(weights, self.num_species, positive=True)
        scalar = _np.ndim(temperature) == 0
        ts = (
            _np.full(self.num_species, float(temperature))
            if scalar
            else _vector(temperature, self.num_species, "temperature")
        )
        if not _np.all(_np.isfinite(ts)) or _np.any(ts <= 0):
            raise ValueError("Temperatures must be positive and finite")
        n = w / w.sum() / self._mw
        if prop_type == ENTHALPY:
            return self._helper.reactant_enthalpy(w.tolist(), ts.tolist())
        if prop_type in (VOLUME, DENSITY, ENTROPY, GIBBS_ENERGY):
            if pressure is None or not _np.isfinite(pressure) or pressure <= 0:
                raise ValueError("A positive pressure in bar is required")
        if prop_type in (VOLUME, DENSITY):
            volume = float(R * _np.dot(n, ts) / (pressure * 100))
            return volume if prop_type == VOLUME else 1 / volume
        if prop_type not in (ENERGY, FROZEN_CP, FROZEN_CV, ENTROPY, GIBBS_ENERGY):
            raise NotImplementedError(f"Mixture.calc_property({prop_type})")
        props = [
            r.thermo(float(t)) if amount > 0 else (0, 0, 0)
            for r, t, amount in zip(self._reactants, ts, w)
        ]
        if prop_type == ENERGY:
            h = self._helper.reactant_enthalpy(w.tolist(), ts.tolist())
            return float(h - R * _np.dot(n[self._gas], ts[self._gas]))
        if prop_type in (FROZEN_CP, FROZEN_CV):
            if any(p[0] is None for p in props):
                raise NotImplementedError(
                    "Mixture heat capacity for reference-only reactants"
                )
            cp = float(_np.dot(n, [p[0] for p in props]) * 1000)
            return cp if prop_type == FROZEN_CP else cp - float(n[self._gas].sum() * R)
        if any(p[2] is None for p in props):
            raise NotImplementedError("Mixture entropy for reference-only reactants")
        ng = n[self._gas].sum()
        entropy_terms = _np.array([p[2] * 1000 for p in props])
        active = self._gas & (n > 1e-35)
        entropy_terms[active] -= R * _np.log(n[active] * pressure / ng)
        entropy = float(_np.dot(n, entropy_terms))
        if prop_type == ENTROPY:
            return entropy if scalar else entropy / 1000
        h = _np.array([p[1] * 1000 for p in props])
        g = float(_np.dot(n, h - ts * entropy_terms))
        # NASA's array branch uses P/1e5 with bar input and unnormalized masses.
        return (
            g
            if scalar
            else (g - R * float(_np.dot(n[self._gas], ts[self._gas])) * _np.log(1e5))
            * float(w.sum())
        )


class _Solver:
    def __init__(
        self,
        products,
        *,
        reactants=None,
        transport=False,
        ions=False,
        trace=-1.0,
        insert=[],
        smooth_truncation=False,
        truncation_width=-1.0,
    ):
        if not isinstance(products, Mixture) or (
            reactants is not None and not isinstance(reactants, Mixture)
        ):
            raise TypeError("products and reactants must be Mixture objects")
        if smooth_truncation or truncation_width > 0:
            raise NotImplementedError("smooth_truncation / truncation_width")
        self._products = products
        self._reactants = products if reactants is None else reactants
        self._transport = bool(transport)
        self._inner = _rs.Solver(
            self._reactants._reactants,
            ions=ions,
            only=products.species_names,
            insert=list(insert),
            trace=0 if trace < 0 else trace,
        )

    @property
    def num_reactants(self):
        return self._reactants.num_species

    @property
    def num_products(self):
        return self._products.num_species

    @property
    def num_gas(self):
        return int(self._products._gas.sum())

    @property
    def num_condensed(self):
        return self.num_products - self.num_gas

    @property
    def num_elements(self):
        return len(set().union(*(r._nasa_formula for r in self._products._reactants)))

    def _begin(self, soln):
        if soln._solver is not self:
            raise ValueError("Solution belongs to a different solver")
        soln._states = []
        soln._extra = {}
        soln.last_error = NOT_CONVERGED

    def _finish(self, soln, states, transports=None):
        values = (
            transports
            if transports is not None
            else [s.transport() if self._transport else None for s in states]
        )
        soln._transports = values
        soln._states = states
        soln.last_error = SUCCESS


class EqSolver(_Solver):
    @property
    def max_equations(self):
        return self.num_elements + self.num_condensed + 2

    def solve(self, soln, eq_type, state1, state2, amounts):
        self._begin(soln)
        w = _vector(amounts, self.num_reactants, positive=True)
        state = self._inner.solve(*_state_inputs(eq_type, state1, state2), w.tolist())
        self._finish(soln, [state])
        soln._inputs = (eq_type, float(state1), float(state2), w.copy())


class RocketSolver(_Solver):
    def solve(
        self,
        soln,
        weights,
        pc,
        pi_p=None,
        subar=None,
        supar=None,
        iac=True,
        n_frz=None,
        hc=None,
        tc=None,
        mdot=None,
        ac_at=None,
        tc_est=None,
    ):
        self._begin(soln)
        w = _vector(weights, self.num_reactants, positive=True)
        if (hc is None) == (tc is None):
            raise ValueError("RocketSolver.solve requires exactly one of hc or tc")
        if n_frz is not None and n_frz < 0:
            raise NotImplementedError(
                "RocketSolver n_frz=-1: frozen reactant composition"
            )
        for name, value in (("pi_p", pi_p), ("subar", subar), ("supar", supar)):
            if value is not None and _np.ndim(value) != 0 and len(value) == 0:
                raise ValueError(f"{name} cannot be empty")
        result = self._inner.rocket(
            w.tolist(),
            [298.15] * len(w),
            _scalar32(pc) * 1e5,
            pi_p=pi_p,
            subar=subar,
            supar=supar,
            iac=iac,
            n_frz=n_frz,
            hc=None if hc is None else float(hc) * R,
            tc=tc,
            mdot=mdot,
            ac_at=ac_at,
            tc_est=tc_est,
            transport=self._transport,
        )
        stations = result.stations
        chamber = None
        if iac and n_frz == 1:
            chamber = self._inner.solve(
                "TP" if tc is not None else "HP",
                float(tc) if tc is not None else float(hc) * R,
                _scalar32(pc) * 1e5,
                w.tolist(),
            )
        self._finish(soln, [s.state for s in stations], [s.transport for s in stations])
        if chamber is not None:
            values = (
                _np.array(
                    [
                        s.state.cp_equilibrium_j_kg_k
                        if s.state.equilibrium
                        else s.state.cp_frozen_j_kg_k
                        for s in stations
                    ]
                )
                / 1000
            )
            values[0] = chamber.cp_equilibrium_j_kg_k / 1000
            soln._extra["cp_eq"] = values
        for name, attr in {
            "Mach": "mach",
            "sonic_velocity": "sound_speed_m_s",
            "ae_at": "area_ratio",
            "c_star": "cstar_m_s",
            "coefficient_of_thrust": "thrust_coefficient",
            "Isp": "isp_m_s",
            "Isp_vacuum": "isp_vacuum_m_s",
        }.items():
            soln._extra[name] = _np.array([getattr(s, attr) or 0.0 for s in stations])
        if not result.complete:
            soln.last_error = NOT_CONVERGED
            _warnings.warn(
                f"RocketSolver.solve: {result.stop_reason}",
                RuntimeWarning,
                stacklevel=2,
            )


class ShockSolver(_Solver):
    def solve(
        self,
        soln,
        weights,
        T0,
        p0,
        u1=None,
        Mach1=None,
        reflected=True,
        incident_frozen=False,
        reflected_frozen=False,
    ):
        self._begin(soln)
        w = _vector(weights, self.num_reactants, positive=True)
        if bool(reflected) != soln._reflected:
            raise ValueError("ShockSolution reflected must match ShockSolver.solve")
        result = self._inner.shock(
            w.tolist(),
            _scalar32(T0),
            _scalar32(p0) * 1e5,
            u1=u1,
            mach1=Mach1,
            reflected=reflected,
            incident_frozen=incident_frozen,
            reflected_frozen=reflected_frozen,
        )
        stations = [result.station1, result.station2]
        if result.station5 is not None:
            stations.append(result.station5)
        states = [s.state for s in stations]
        cv = _np.array(
            [s.cv_equilibrium_j_kg_k / 1000 if s.equilibrium else 0.0 for s in states]
        )
        self._finish(
            soln,
            states,
            [None]
            + [s.state.transport() if self._transport else None for s in stations[1:]],
        )
        soln._extra["cv_eq"] = cv
        for name in ("conductivity_eq", "Pr_eq"):
            attr, scale = _TRANSPORT_FIELDS[name]
            soln._extra[name] = _np.array(
                [
                    getattr(tr, attr) * scale
                    if tr is not None and state.equilibrium
                    else 0.0
                    for state, tr in zip(soln._states, soln._transports)
                ]
            )
        for name, attr in {
            "velocity": "velocity_m_s",
            "Mach": "mach",
            "sonic_velocity": "sonic_velocity_m_s",
        }.items():
            soln._extra[name] = _np.array([getattr(s, attr) for s in stations])
        soln._extra.update(
            rho12=1 / result.rho2_rho1,
            rho52=result.rho5_rho2 or 0.0,
            P21=result.p2_p1,
            P52=result.p5_p2 or 0.0,
            T21=result.t2_t1,
            T52=result.t5_t2 or 0.0,
            M21=result.m2_m1,
            M52=result.m5_m2 or 0.0,
            v2=result.v2_m_s,
            u5_p_v2=result.u5_plus_v2_m_s or 0.0,
        )


class DetonationSolver(_Solver):
    def solve(self, soln, weights, T1, p1, frozen=False):
        self._begin(soln)
        w = _vector(weights, self.num_reactants, positive=True)
        result = self._inner.detonation(
            w.tolist(), _scalar32(T1), _scalar32(p1) * 1e5, frozen=frozen
        )
        self._finish(soln, [result.burned])
        u = result.unburned
        soln._extra.update(
            P1=u.pressure_pa / 1e5,
            T1=u.temperature_k,
            H1=u.enthalpy_j_kg / 1000,
            M1=u.molecular_weight_kg_kmol,
            gamma1=u.gamma_s,
            sonic_velocity1=u.sound_speed_m_s,
            Mach=result.mach,
            velocity=result.velocity_m_s,
            sonic_velocity=result.burned.sound_speed_m_s,
            P_P1=result.pressure_ratio,
            T_T1=result.temperature_ratio,
            M_M1=result.molecular_weight_ratio,
            rho_rho1=result.density_ratio,
        )


_STATE_FIELDS = {
    "T": ("temperature_k", 1),
    "P": ("pressure_pa", 1e-5),
    "volume": ("specific_volume_m3_kg", 1),
    "density": ("density_kg_m3", 1),
    "M": ("M", 1),
    "MW": ("MW", 1),
    "enthalpy": ("enthalpy_j_kg", 0.001),
    "energy": ("energy_j_kg", 0.001),
    "entropy": ("entropy_j_kg_k", 0.001),
    "gibbs_energy": ("gibbs_energy_j_kg", 0.001),
    "gamma_s": ("gamma_s", 1),
    "cp_fr": ("cp_frozen_j_kg_k", 0.001),
    "cp_eq": ("cp_equilibrium_j_kg_k", 0.001),
    "cv_fr": ("cv_frozen_j_kg_k", 0.001),
    "cv_eq": ("cv_equilibrium_j_kg_k", 0.001),
}
_TRANSPORT_FIELDS = {
    "viscosity": ("viscosity_pa_s", 1e4),
    "conductivity_fr": ("conductivity_frozen_w_m_k", 10),
    "conductivity_eq": ("conductivity_equilibrium_w_m_k", 10),
    "Pr_fr": ("prandtl_frozen", 1),
    "Pr_eq": ("prandtl_equilibrium", 1),
}


def _field(name):
    return property(lambda self: self._get(name))


class _Solution:
    _array = False

    def __init__(self, solver):
        self._solver = solver
        self._states = []
        self._transports = []
        self._extra = {}
        self.last_error = SUCCESS

    @property
    def converged(self):
        return bool(self._states) and self.last_error == SUCCESS

    def _shape(self, values):
        return _np.array(values, dtype=_np.float64) if self._array else float(values[0])

    def _get(self, name):
        if not self._states:
            raise RuntimeError(f"{type(self).__name__}.{name}: no verified solution")
        name = {"cp": "cp_eq", "cv": "cv_eq"}.get(name, name)
        if name in self._extra:
            value = self._extra[name]
            return value.copy() if isinstance(value, _np.ndarray) else value
        if name in _TRANSPORT_FIELDS:
            attr, scale = _TRANSPORT_FIELDS[name]
            return self._shape(
                [
                    0.0 if t is None else getattr(t, attr) * scale
                    for t in self._transports
                ]
            )
        attr, scale = _STATE_FIELDS[name]
        values = []
        for state, tr in zip(self._states, self._transports):
            if name == "cp_fr" and tr is not None:
                value = tr.cp_frozen_j_kg_k
            elif name == "cp_eq" and not state.equilibrium:
                value = state.cp_frozen_j_kg_k
            elif name == "cv_eq" and not state.equilibrium:
                value = state.cv_frozen_j_kg_k
            elif name == "cv_eq" and state.cp_equilibrium_j_kg_k == 0:
                value = 0.0
            else:
                value = getattr(state, attr)
            values.append(value * scale)
        return self._shape(values)

    @property
    def n(self):
        return (
            self._shape([1 / s.M for s in self._states])
            if self._states
            else self._get("n")
        )

    @property
    def nj(self):
        if not self._states:
            raise RuntimeError("nj: no verified solution")
        rows = [
            [s._nasa_nj.get(n, 0.0) for n in self._solver._products._names]
            for s in self._states
        ]
        return _np.array(rows if self._array else rows[0], dtype=_np.float64)

    @property
    def ln_nj(self):
        if not self._states:
            raise RuntimeError("ln_nj: no verified solution")
        logs = [s._nasa_ln_nj for s in self._states]
        if any(row is None for row in logs):
            raise NotImplementedError("untruncated ln_nj for frozen stations")
        names = [
            n
            for n, gas in zip(
                self._solver._products._names, self._solver._products._gas
            )
            if gas
        ]
        if any(n not in row for row in logs for n in names):
            raise NotImplementedError("ln_nj for products removed by element reduction")
        rows = [[row[n] for n in names] for row in logs]
        return _np.array(rows if self._array else rows[0], dtype=_np.float64)

    def _fractions(self, name):
        if not self._states:
            raise RuntimeError(f"{name}: no verified solution")
        rows = [getattr(s, name) for s in self._states]
        return {
            n: self._shape([row.get(n, 0.0) for row in rows])
            for n in self._solver._products._names
        }

    @property
    def mass_fractions(self):
        return self._fractions("mass_fractions")

    @property
    def mole_fractions(self):
        return self._fractions("mole_fractions")


class EqSolution(_Solution):
    def __init__(self, solver, T_init=None, nj_init=None):
        if T_init is not None or nj_init is not None:
            raise NotImplementedError("EqSolution T_init / nj_init initial guesses")
        super().__init__(solver)


class RocketSolution(_Solution):
    _array = True

    @property
    def num_pts(self):
        return len(self._states)


class ShockSolution(_Solution):
    _array = True

    def __init__(self, solver, reflected=True):
        super().__init__(solver)
        self._reflected = bool(reflected)


class DetonationSolution(_Solution):
    pass


for _cls in (EqSolution, RocketSolution, ShockSolution, DetonationSolution):
    for _name in (*_STATE_FIELDS, *_TRANSPORT_FIELDS):
        if _cls is DetonationSolution and _name == "volume":
            continue
        setattr(_cls, _name, _field(_name))
    if _cls is not DetonationSolution:
        for _name in ("cp", "cv"):
            setattr(_cls, _name, _field(_name))
for _cls, _names in (
    (
        RocketSolution,
        (
            "Mach",
            "sonic_velocity",
            "ae_at",
            "c_star",
            "coefficient_of_thrust",
            "Isp",
            "Isp_vacuum",
        ),
    ),
    (
        ShockSolution,
        (
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
    ),
    (
        DetonationSolution,
        (
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
    ),
):
    for _name in _names:
        setattr(_cls, _name, _field(_name))


class EqDerivatives:
    """NASA-shaped sensitivities; state inputs follow EqSolver, outputs use kJ."""

    def __init__(self, solver, solution):
        if solution._solver is not solver:
            raise ValueError("Solution belongs to a different solver")
        self._solver = solver
        self._solution = solution
        self._results = {}
        self.last_error = SUCCESS

    def compute_derivatives(self, check_closure_defect=False):
        if not self._solution.converged:
            raise ValueError("EqDerivatives needs a converged solution")
        eq_type, first, second, w = self._solution._inputs
        mode, a, b = _state_inputs(eq_type, first, second)
        self._results.clear()
        self.last_error = NOT_CONVERGED
        d = self._solver._inner.derivatives(mode, a, b, w.tolist())
        scale = _np.array(
            [
                1 if mode.startswith("T") else R,
                1e5 if mode.endswith("P") else 1,
                *([1] * len(w)),
            ]
        )
        rows = {
            "T": _np.array(d.temperature) * scale,
            "n": _np.array(d.total_moles) * scale,
        }
        for symbol, attr in (
            ("H", "enthalpy"),
            ("U", "energy"),
            ("G", "gibbs_energy"),
            ("S", "entropy"),
        ):
            rows[symbol] = _np.array(getattr(d, attr)) * scale / 1000
        rows["nj"] = (
            _np.array(
                [
                    d.species_amounts.get(n, [0.0] * len(scale))
                    for n in self._solver._products._names
                ]
            )
            * scale
        )
        self._store(rows)
        self.last_error = SUCCESS

    def _store(self, rows, suffix=""):
        for symbol, row in rows.items():
            self._results[f"d{symbol}_dstate1{suffix}"] = (
                row[..., 0].copy() if symbol == "nj" else float(row[0])
            )
            self._results[f"d{symbol}_dstate2{suffix}"] = (
                row[..., 1].copy() if symbol == "nj" else float(row[1])
            )
            self._results[f"d{symbol}_dw0{suffix}"] = row[..., 2:].copy()

    def compute_fd(self, h, verbose=False, central=False):
        if verbose:
            raise NotImplementedError("EqDerivatives.compute_fd verbose output")
        if not self._solution.converged or not _np.isfinite(h) or h <= 0:
            raise ValueError(
                "compute_fd needs a converged solution and positive finite h"
            )
        eq_type, a, b, w = self._solution._inputs
        inputs = _np.array([a, b, *w])
        rows = {s: [] for s in ("T", "n", "H", "U", "G", "S", "nj")}

        def sample(v):
            x = EqSolution(self._solver)
            self._solver.solve(x, eq_type, v[0], v[1], v[2:])
            return {
                "T": x.T,
                "n": x.n,
                "H": x.enthalpy,
                "U": x.energy,
                "G": x.gibbs_energy,
                "S": x.entropy,
                "nj": x.nj,
            }

        self.last_error = NOT_CONVERGED
        base = sample(inputs)
        for i, value in enumerate(inputs):
            step = float(h) * max(1, abs(value))
            hi = inputs.copy()
            hi[i] += step
            plus = sample(hi)
            lo = inputs.copy()
            lo[i] -= step
            minus = sample(lo) if central else base
            for symbol in rows:
                rows[symbol].append(
                    (plus[symbol] - minus[symbol]) / (step * (2 if central else 1))
                )
        self._store({s: _np.array(v).T for s, v in rows.items()}, "_fd")
        self.last_error = SUCCESS

    def _get(self, name):
        if name not in self._results:
            raise RuntimeError(f"EqDerivatives.{name}: compute derivatives first")
        v = self._results[name]
        return v.copy() if isinstance(v, _np.ndarray) else v


for _symbol in ("T", "n", "H", "U", "G", "S", "nj"):
    for _input in ("state1", "state2", "w0"):
        for _suffix in ("", "_fd"):
            _name = f"d{_symbol}_d{_input}{_suffix}"
            setattr(EqDerivatives, _name, _field(_name))


def eq_solve(
    eq_type,
    reactants,
    T=None,
    H=None,
    S=None,
    U=None,
    P=None,
    V=None,
    T_reac=None,
    fuel_amounts=None,
    oxid_amounts=None,
    moles=False,
    of_ratio=None,
    phi=None,
    r_eq=None,
    pct_fuel=None,
    only=None,
    omit=None,
    insert=None,
    trace=None,
    transport=False,
    ions=False,
):
    """Legacy root helper, including 3.3.4's direct H/S/U-to-state1 convention."""
    _warnings.warn(
        "cea.eq_solve is deprecated; use solver objects instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    reac = Mixture(reactants, ions=ions)
    prod = (
        Mixture(only, ions=ions)
        if only
        else Mixture(reactants, True, omit=omit or [], ions=ions)
    )
    if fuel_amounts is None or oxid_amounts is None:
        raise TypeError("Fuel and oxidizer amounts not defined")
    fuel = _vector(fuel_amounts, reac.num_species)
    oxid = _vector(oxid_amounts, reac.num_species)
    if moles:
        fuel, oxid = reac.moles_to_weights(fuel), reac.moles_to_weights(oxid)
    if of_ratio is not None:
        ratio = of_ratio
    elif phi is not None:
        ratio = reac.weight_eq_ratio_to_of_ratio(oxid, fuel, phi)
    elif r_eq is not None:
        ratio = reac.chem_eq_ratio_to_of_ratio(oxid, fuel, r_eq)
    elif pct_fuel is not None:
        ratio = (100 - pct_fuel) / pct_fuel
    else:
        ratio = None
    weights = (
        fuel + oxid if ratio is None else reac.of_ratio_to_weights(oxid, fuel, ratio)
    )
    first = next((v for v in (T, H, S, U) if v is not None), None)
    if first is None:
        if T_reac is None:
            raise TypeError("Reactant temperature not defined")
        if eq_type in (TP, TV):
            first = float(T_reac) if _np.ndim(T_reac) == 0 else T_reac[0]
            _warnings.warn(
                "Problem temperature not defined; using first reactant temperature",
                stacklevel=2,
            )
        else:
            first = reac.calc_property(
                {HP: ENTHALPY, SP: ENTROPY, SV: ENTROPY, UV: ENERGY}[eq_type],
                weights,
                T_reac,
            )
    second = P if P is not None else V
    if second is None:
        raise TypeError("Pressure or volume state not defined")
    solver = EqSolver(
        prod,
        reactants=reac,
        transport=transport,
        ions=ions,
        trace=-1.0 if trace is None else trace,
        insert=insert or [],
    )
    soln = EqSolution(solver)
    solver.solve(soln, eq_type, first, second, weights)
    return soln


__all__ = [
    *constants.__all__,
    "__version__",
    "units",
    "lib_version",
    "lib_version_major",
    "lib_version_minor",
    "lib_version_patch",
    "init",
    "init_thermo",
    "init_trans",
    "is_initialized",
    "set_log_level",
    "Reactant",
    "Mixture",
    "EqSolver",
    "EqSolution",
    "EqDerivatives",
    "RocketSolver",
    "RocketSolution",
    "ShockSolver",
    "ShockSolution",
    "DetonationSolver",
    "DetonationSolution",
    "eq_solve",
]

from . import lib as lib
