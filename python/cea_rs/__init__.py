"""Expose the experimental Rust solver, using SI units throughout."""

from ._cea_rs import (
    ConvergenceError,
    DomainError,
    DetonationResult,
    R,
    Reactant,
    RocketResult,
    ShockResult,
    ShockStation,
    State,
    Station,
    UnburnedState,
    __version__,
)
from ._cea_rs import (
    Solver as _Solver,
)

NASA_THERMO_SHA256 = "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363"

__all__ = [
    "NASA_THERMO_SHA256",
    "ConvergenceError",
    "DomainError",
    "DetonationResult",
    "R",
    "Reactant",
    "RocketResult",
    "ShockResult",
    "ShockStation",
    "Solver",
    "State",
    "Station",
    "UnburnedState",
    "__version__",
    "capabilities",
]


class Solver:
    """Construct an independent Rust solver from named or custom reactants.

    Args:
        reactants: NASA names or Reactant objects, in weight-vector order.
        condensed: Allow condensed equilibrium products.
        ions: Include charged species and electron balance.
        only: Optional explicit product names. This changes the physical model.
        insert: Condensed product names to activate at the initial iteration.
        omit: Product names to exclude from equilibrium.
        trace: NASA trace convergence level; zero selects the default.
    """

    def __init__(
        self,
        reactants: list[str | Reactant],
        *,
        condensed: bool = True,
        ions: bool = False,
        only: list[str] | None = None,
        insert: list[str] | None = None,
        omit: list[str] | None = None,
        trace: float = 0.0,
    ) -> None:
        """Validate the species model before any numerical solve."""
        self._inner = _Solver(
            [Reactant(r) if isinstance(r, str) else r for r in reactants],
            condensed=condensed,
            ions=ions,
            only=only,
            insert=insert or [],
            omit=omit or [],
            trace=trace,
        )

    def solve(
        self, mode: str, state1: float, state2: float, weights: list[float]
    ) -> State:
        """Solve TP, HP, SP, TV, UV or SV with SI state inputs.

        Args:
            mode: One of TP, HP, SP, TV, UV, SV (case-sensitive).
            state1: Temperature K, enthalpy J/kg, entropy J/kg/K, or energy J/kg.
            state2: Pressure Pa for TP/HP/SP, specific volume m³/kg otherwise.
            weights: Nonnegative reactant masses, normalized internally.

        Returns:
            A verified state. Failed solves raise ConvergenceError; fitted-data
            domain violations raise DomainError. No stale state is returned.
            Mole fractions follow NASA's reporting convention, with condensed
            moles included in the denominator; they are not gas-only fractions.
            Frozen states expose equilibrium Cp as None. Out-of-fit trace
            species at <= 1e-11 mole fraction are listed explicitly.
        """
        return self._inner.solve(mode, state1, state2, weights)

    def reactant_enthalpy(
        self, weights: list[float], temperatures_k: list[float]
    ) -> float:
        """Compute stream enthalpy in J/kg at each reactant's temperature."""
        return self._inner.reactant_enthalpy(weights, temperatures_k)

    def rocket(
        self,
        weights: list[float],
        temperatures_k: list[float],
        chamber_pressure_pa: float,
        expansion_ratio: float | None = None,
        *,
        flow: str = "equilibrium",
        iac: bool = True,
        transport: bool = False,
        pi_p: float | list[float] | None = None,
        subar: float | list[float] | None = None,
        supar: float | list[float] | None = None,
        ac_at: float | None = None,
        mdot: float | None = None,
        n_frz: int | None = None,
        hc: float | None = None,
        tc: float | None = None,
        tc_est: float | None = None,
    ) -> RocketResult:
        """Solve chamber, throat and pressure/area exit schedules in SI.

        The existing positional expansion_ratio specifies one supersonic exit.
        Otherwise pi_p, subar and supar accept scalars or lists, processed in
        that order. pi_p is injector/chamber pressure divided by exit pressure.
        FAC (iac=False) requires ac_at or mdot in kg/(s m²), and returns injector,
        infinity, combustor end and throat before the requested exits.

        n_frz counts stations from one, including FAC's infinity station; zero
        means equilibrium. Points upstream of the freeze point are omitted as
        in NASA CEA. A frozen fit-limit stop returns complete=False and only
        verified stations, with stop_reason explaining the limit.

        hc (J/kg) or tc (K) assigns the stagnation state instead of the reactant
        enthalpy. tc_est seeds a chamber temperature iteration. transport must be
        False. Specific impulses are velocities in m/s; divide by 9.80665 to obtain seconds.
        """

        def schedule(value: float | list[float] | None) -> list[float]:
            if value is None:
                return []
            if isinstance(value, (float, int)):
                return [float(value)]
            return list(value)

        return self._inner.rocket(
            weights,
            temperatures_k,
            chamber_pressure_pa,
            expansion_ratio,
            flow=flow,
            iac=iac,
            transport=transport,
            pi_p=schedule(pi_p),
            subar=schedule(subar),
            supar=schedule(supar),
            ac_at=ac_at,
            mdot=mdot,
            n_frz=n_frz,
            hc=hc,
            tc=tc,
            tc_est=tc_est,
        )

    def shock(
        self,
        weights: list[float],
        temperature_k: float,
        pressure_pa: float,
        *,
        u1: float | None = None,
        mach1: float | None = None,
        reflected: bool = False,
        incident_frozen: bool = False,
        reflected_frozen: bool = False,
    ) -> ShockResult:
        """Solve a NASA gas-mixture shock in SI, using exactly one of u1/mach1.

        Station velocities use each shock's frame; v2_m_s uses the wall frame.
        Any incident/reflected frozen/equilibrium combination is supported.
        Reflected-frozen sound speed and gamma retain NASA's station-2 Cp basis;
        cp_frozen_j_kg_k is evaluated at the actual station temperature.
        Failed solves raise an error and never return a retained earlier state.
        """
        return self._inner.shock(
            weights, temperature_k, pressure_pa, u1=u1, mach1=mach1,
            reflected=reflected, incident_frozen=incident_frozen,
            reflected_frozen=reflected_frozen,
        )

    def detonation(
        self,
        weights: list[float],
        temperature_k: float,
        pressure_pa: float,
        *,
        frozen: bool = False,
    ) -> DetonationResult:
        """Solve a Chapman-Jouguet wave from unburned gas temperature and pressure.

        Inputs and outputs use SI; weights are reactant masses. The unburned
        composition stays fixed and the burned state is in equilibrium.
        Reactants must be named gas species with NASA fits. frozen=True raises
        NotImplementedError here; native CEA 3.3.4 also leaves it unimplemented.
        """
        return self._inner.detonation(
            weights, temperature_k, pressure_pa, frozen=frozen
        )


def capabilities() -> dict[str, object]:
    """Describe implemented algorithms separately from established parity."""
    return {
        "equilibrium": ["TP", "HP", "SP", "TV", "UV", "SV"],
        "rocket": [
            "IAC equilibrium",
            "IAC frozen at chamber",
            "IAC frozen at throat",
            "FAC",
            "pressure and area station schedules",
            "station-index freezing",
        ],
        "condensed": "experimental active-set equilibrium",
        "ions": "experimental charge-neutral equilibrium",
        "unsupported": [
            "shock",
            "detonation",
            "transport",
            "legacy input decks",
            "NASA Python API compatibility",
            "analytic sensitivities",
        ],
        "full_cea_parity": False,
        "production_default": "configured by consumer",
    }
