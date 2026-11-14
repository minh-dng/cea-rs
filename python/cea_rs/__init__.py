"""Expose the experimental Rust solver, using SI units throughout."""

from ._cea_rs import (
    ConvergenceError,
    DomainError,
    EqDerivatives,
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
    "EqDerivatives",
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

    def of_ratio_to_weights(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> list[float]:
        """Return NASA's scaled reactant masses for an oxidant/fuel mass ratio.

        Each stream vector follows the full reactant order. Its internal ratios
        are normalized separately; the returned vector need not sum to one.
        For percent fuel p, use of_ratio=(100-p)/p, with 0 < p <= 100.
        """
        return self._inner.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

    def chem_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: list[float],
        fuel_weights: list[float],
        chem_eq_ratio: float,
    ) -> float:
        """Convert chemical equivalence ratio r using NASA's element valences."""
        return self._inner.chem_eq_ratio_to_of_ratio(
            oxidant_weights, fuel_weights, chem_eq_ratio
        )

    def of_ratio_to_chem_eq_ratio(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> float:
        return self._inner.of_ratio_to_chem_eq_ratio(
            oxidant_weights, fuel_weights, of_ratio
        )

    def weight_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: list[float],
        fuel_weights: list[float],
        weight_eq_ratio: float,
    ) -> float:
        """Convert mass equivalence ratio phi to oxidant/fuel mass ratio."""
        return self._inner.weight_eq_ratio_to_of_ratio(
            oxidant_weights, fuel_weights, weight_eq_ratio
        )

    def of_ratio_to_weight_eq_ratio(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> float:
        return self._inner.of_ratio_to_weight_eq_ratio(
            oxidant_weights, fuel_weights, of_ratio
        )

    def weights_to_moles(self, weights: list[float]) -> list[float]:
        """Convert kg to kmol without normalizing."""
        return self._inner.weights_to_moles(weights)

    def moles_to_weights(self, moles: list[float]) -> list[float]:
        """Convert kmol to kg without normalizing."""
        return self._inner.moles_to_weights(moles)

    def get_valence(self) -> dict[str, tuple[float, float]]:
        """Return each element's (negative, positive) NASA valences."""
        return self._inner.get_valence()

    def derivatives(
        self,
        mode: str,
        state1: float,
        state2: float,
        weights: list[float],
        *,
        finite_difference: bool = False,
        relative_step: float = 1e-4,
    ) -> EqDerivatives:
        """Solve afresh and return response derivatives in SI.

        Every response vector is ordered [state1, state2, weights...]; species
        amounts and total gas moles use kmol/kg. Weights are differentiated as
        supplied, including their normalization, while state inputs stay fixed.
        Analytic derivatives assume default hard truncation and fixed phases;
        ion sensitivities require finite_difference=True. Smooth truncation
        and recovery solves with a changed cutoff need additional core data.
        Analytic UV/SV with condensed products is unavailable; finite differences
        remain available wherever the public forward solver converges.
        Finite differences use a central stencil (second-order forward at a
        zero weight) and reject changes to the product pool or active phases.
        """
        return self._inner.derivatives(
            mode,
            state1,
            state2,
            weights,
            finite_difference=finite_difference,
            relative_step=relative_step,
        )

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
