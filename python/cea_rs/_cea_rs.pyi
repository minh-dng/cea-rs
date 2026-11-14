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
    def thermo(
        self, temperature_k: float
    ) -> tuple[float | None, float, float | None]: ...

class State:
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
    name: str
    sound_speed_m_s: float
    cstar_m_s: float
    thrust_coefficient: float
    isp_m_s: float
    isp_vacuum_m_s: float
    state: State
    velocity_m_s: float
    mach: float
    area_ratio: float | None

class RocketResult:
    stations: list[Station]
    throat_index: int
    complete: bool
    stop_reason: str | None
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
        pi_p: list[float] = ...,
        subar: list[float] = ...,
        supar: list[float] = ...,
        ac_at: float | None = None,
        mdot: float | None = None,
        n_frz: int | None = None,
        hc: float | None = None,
        tc: float | None = None,
        tc_est: float | None = None,
    ) -> RocketResult: ...
