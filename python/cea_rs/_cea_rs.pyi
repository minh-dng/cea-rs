"""Type the Rust extension's SI interface."""

class ConvergenceError(RuntimeError): ...
class DomainError(ValueError): ...

R: float
__version__: str

class Reactant:
    name: str
    molecular_weight: float
    def __init__(
        self,
        name: str,
        *,
        formula: dict[str, float] | None = None,
        molecular_weight: float | None = None,
        enthalpy_j_mol: float | None = None,
        reference_temperature_k: float = 298.15,
    ) -> None: ...
    def get_valid_temperature_range(self) -> tuple[float, float]: ...
    def thermo(
        self, temperature_k: float
    ) -> tuple[float | None, float, float | None]: ...

class State:
    gibbs_energy_j_kg: float
    cv_frozen_j_kg_k: float
    cv_equilibrium_j_kg_k: float | None
    specific_volume_m3_kg: float
    M: float
    MW: float
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    enthalpy_j_kg: float
    energy_j_kg: float
    entropy_j_kg_k: float
    cp_equilibrium_j_kg_k: float | None
    cp_frozen_j_kg_k: float
    gamma_s: float
    molecular_weight_kg_kmol: float
    sound_speed_m_s: float
    equilibrium: bool
    iterations: int
    element_residual: float
    constraint_residual: float
    trace_extrapolated_species: list[str]
    mole_fractions: dict[str, float]
    mass_fractions: dict[str, float]

class Station:
    state: State
    velocity_m_s: float
    mach: float
    area_ratio: float | None

class RocketResult:
    chamber: Station
    throat: Station
    exit: Station
    cstar_m_s: float
    vacuum_velocity_m_s: float

class Solver:
    def __init__(
        self,
        reactants: list[Reactant],
        *,
        condensed: bool = True,
        ions: bool = False,
        only: list[str] | None = None,
        insert: list[str] = ...,
    ) -> None: ...
    def solve(
        self, mode: str, state1: float, state2: float, weights: list[float]
    ) -> State: ...
    def reactant_enthalpy(
        self, weights: list[float], temperatures_k: list[float]
    ) -> float: ...
    def of_ratio_to_weights(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> list[float]: ...
    def chem_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: list[float],
        fuel_weights: list[float],
        chem_eq_ratio: float,
    ) -> float: ...
    def of_ratio_to_chem_eq_ratio(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> float: ...
    def weight_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: list[float],
        fuel_weights: list[float],
        weight_eq_ratio: float,
    ) -> float: ...
    def of_ratio_to_weight_eq_ratio(
        self, oxidant_weights: list[float], fuel_weights: list[float], of_ratio: float
    ) -> float: ...
    def weights_to_moles(self, weights: list[float]) -> list[float]: ...
    def moles_to_weights(self, moles: list[float]) -> list[float]: ...
    def get_valence(self) -> dict[str, tuple[float, float]]: ...
    def derivatives(
        self,
        mode: str,
        state1: float,
        state2: float,
        weights: list[float],
        *,
        finite_difference: bool = False,
        relative_step: float = 1e-4,
    ) -> EqDerivatives: ...
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
    ) -> RocketResult: ...

class EqDerivatives:
    """SI response vectors in order [state1, state2, reactant weights...]."""
    temperature: list[float]
    total_moles: list[float]
    species_amounts: dict[str, list[float]]
    enthalpy: list[float]
    energy: list[float]
    gibbs_energy: list[float]
    entropy: list[float]
    closure_defect: float | None
