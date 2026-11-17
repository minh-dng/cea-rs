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
    extrapolated_species: list[str]
    mole_fractions: dict[str, float]
    mass_fractions: dict[str, float]
    def transport(self) -> Transport: ...

class Transport:
    viscosity_pa_s: float
    conductivity_frozen_w_m_k: float
    conductivity_equilibrium_w_m_k: float
    prandtl_frozen: float
    prandtl_equilibrium: float

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
    ) -> ShockResult: ...
    def __init__(
        self,
        reactants: list[Reactant],
        *,
        condensed: bool = True,
        ions: bool = False,
        only: list[str] | None = None,
        insert: list[str] = ...,
        omit: list[str] = ...,
        trace: float = 0.0,
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
    def detonation(
        temperature_k: float,
        pressure_pa: float,
        frozen: bool = False,
    ) -> DetonationResult: ...
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
class UnburnedState:
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    enthalpy_j_kg: float
    molecular_weight_kg_kmol: float
    gamma_s: float
    sound_speed_m_s: float

class DetonationResult:
    unburned: UnburnedState
    burned: State
    pressure_ratio: float
    temperature_ratio: float
    molecular_weight_ratio: float
    density_ratio: float
    mach: float
    velocity_m_s: float
class ShockStation:
    state: State
    sonic_velocity_m_s: float

class ShockResult:
    station1: ShockStation
    station2: ShockStation
    station5: ShockStation | None
    p2_p1: float
    t2_t1: float
    m2_m1: float
    rho2_rho1: float
    v2_m_s: float
    p5_p2: float | None
    t5_t2: float | None
    m5_m2: float | None
    rho5_rho2: float | None
    u5_plus_v2_m_s: float | None
