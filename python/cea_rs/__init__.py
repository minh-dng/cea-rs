"""Expose the experimental Rust solver, using SI units throughout."""

from ._cea_rs import (
    ConvergenceError,
    DomainError,
    R,
    Reactant,
    RocketResult,
    State,
    Station,
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
    "R",
    "Reactant",
    "RocketResult",
    "Solver",
    "State",
    "Station",
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
    """

    def __init__(
        self,
        reactants: list[str | Reactant],
        *,
        condensed: bool = True,
        ions: bool = False,
        only: list[str] | None = None,
        insert: list[str] | None = None,
    ) -> None:
        """Validate the species model before any numerical solve."""
        self._inner = _Solver(
            [Reactant(r) if isinstance(r, str) else r for r in reactants],
            condensed=condensed,
            ions=ions,
            only=only,
            insert=insert or [],
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
        expansion_ratio: float,
        *,
        flow: str = "equilibrium",
        iac: bool = True,
        transport: bool = False,
    ) -> RocketResult:
        """Solve an infinite-area chamber, sonic throat, and supersonic exit.

        Args:
            weights: Nonnegative reactant masses, normalized internally.
            temperatures_k: Temperatures in reactant order, without extrapolation.
            chamber_pressure_pa: Stagnation chamber pressure in Pa.
            expansion_ratio: Geometric exit area divided by throat area, >= 1.
            flow: equilibrium, frozen_at_chamber, or frozen_at_throat.
            iac: Must be True; finite-area chambers are unimplemented.
            transport: Must be False; transport is unimplemented.

        Returns:
            Three stations and ideal c-star/vacuum velocity in m/s. Divide
            vacuum_velocity_m_s by 9.80665 for vacuum specific impulse in seconds.
        """
        return self._inner.rocket(
            weights,
            temperatures_k,
            chamber_pressure_pa,
            expansion_ratio,
            flow=flow,
            iac=iac,
            transport=transport,
        )


def capabilities() -> dict[str, object]:
    """Describe implemented algorithms separately from established parity."""
    return {
        "equilibrium": ["TP", "HP", "SP", "TV", "UV", "SV"],
        "rocket": ["IAC equilibrium", "IAC frozen at chamber", "IAC frozen at throat"],
        "condensed": "experimental active-set equilibrium",
        "ions": "experimental charge-neutral equilibrium",
        "unsupported": [
            "finite-area chamber",
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
