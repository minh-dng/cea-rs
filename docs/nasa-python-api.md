# NASA Python 3.3.4 API contract

Inspected the PyPI wheel and Python/Cython sources at NASA commit
[2f79a647](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da/source/bind/python)
on 2026-11-17. This inventory records public API facts, not copied implementation.
NASA is Apache-2.0; see [provenance](provenance.md) and the root NOTICE.

## Conventions

- Solver pressure inputs and solution pressures: bar. Temperatures: K.
- EqSolver state1: T for TP/TV, H/R for HP, S/R for SP/SV, U/R for UV,
  where H/U are J/kg, S is J/(kg K), and R=8314.51 J/(kmol K).
  State2 is bar for TP/HP/SP and m³/kg for TV/UV/SV.
- Rocket hc uses H/R despite its docstring saying (kJ/kg)/R. mdot is kg/(s m²).
- Solutions: H/U/G in kJ/kg, Cp/Cv/S in kJ/(kg K), volume m³/kg, density kg/m³,
  M/MW kg/kmol, velocities/c-star/Isp in m/s, viscosity millipoise,
  conductivity µW/(cm K); ratios, fractions and Prandtl numbers are dimensionless.
- Mixture.calc_property: H/U/G J/kg, scalar-temperature S J/(kg K), Cp/Cv J/(kg K).
  In 3.3.4 array-temperature S is divided by 1000; array-temperature G retains
  the supplied weight scale; V/DENSITY retain legacy ×1000/×0.001 relative to SI.
- Weights are masses normalized by the solver; mole inputs require moles_to_weights.
  Conversion helpers preserve amount scale and return float64 NumPy arrays.
- Products are selected from all reactant elements, including zero-weight reactants.
  omit applies only to generated products; explicit Mixture names are the only pool.
  insert seeds condensed phases. ions must be requested in mixtures and solver.
  trace<0 selects the native default; smooth truncation is a distinct unsupported model.
- Eq/detonation thermodynamic properties are float scalars; rocket/shock station
  properties are float64 arrays. nj is (species,) or (stations,species), on a kmol/kg
  basis; ln_nj includes gas species only and retains untruncated concentrations.
  n is total gas kmol/kg. Fraction dictionaries contain every product, including zeros.
- Rocket stations follow chamber/throat/pressure/subsonic/supersonic schedules;
  FAC includes injector/infinity/combustor/throat, with NASA display omissions.
  n_frz is one-based and counts the FAC infinity station; 0 means equilibrium.
  Shock stations are [1,2] or [1,2,5]; ratio properties are scalar floats.
- All solve methods mutate the solution and return None. converged is bool;
  last_error is an integer code. This port raises on failed verification and clears
  the old solution, instead of exposing NASA retained or nonconverged values.

## Classes, methods and properties

Constructor rows expand documented keyword defaults accepted by native kwargs.
Public signatures below omit Cython implementation types.

### Reactant

Constructor: `Reactant(name, formula=None, molecular_weight=None, enthalpy=None, enthalpy_units=None, temperature=None)`.

| Member | Signature / return convention |
|---|---|
| enthalpy | See conventions above |
| enthalpy_units | See conventions above |
| formula | See conventions above |
| get_valid_temperature_range | `get_valid_temperature_range(self)` |
| molecular_weight | See conventions above |
| name | See conventions above |
| temperature | See conventions above |

### Mixture

Constructor: `Mixture(species=None, products_from_reactants=False, omit=[], ions=False)`.

| Member | Signature / return convention |
|---|---|
| calc_property | `calc_property(self, prop_type, weights, temperature: 'float | list | np.ndarray', pressure: 'Optional[float]' = None)` |
| chem_eq_ratio_to_of_ratio | `chem_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, chem_eq_ratio)` |
| moles_to_weights | `moles_to_weights(self, moles)` |
| num_species | int |
| of_ratio_to_weights | `of_ratio_to_weights(self, oxidant_weights, fuel_weights, of_ratio)` |
| species_names | Get list of species names in the mixture. |
| weight_eq_ratio_to_of_ratio | `weight_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, weight_eq_ratio)` |
| weights_to_moles | `weights_to_moles(self, weights)` |

### EqSolver

Constructor: `EqSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`.

| Member | Signature / return convention |
|---|---|
| max_equations | Maximum number of equations in system; num_elements + num_condensed + (1 for TP/TV problems, 2 otherwise). |
| num_condensed | Number of condensed phase species. |
| num_elements | Number of chemical elements. |
| num_gas | Number of gas phase species. |
| num_products | Number of product species. |
| num_reactants | Number of reactant species. |
| solve | `solve(self, soln, eq_type, state1, state2, amounts)` |

### EqSolution

Constructor: `EqSolution(solver, T_init=None, nj_init=None)`.

| Member | Signature / return convention |
|---|---|
| M | Molecular weight, gasses only (1/n). |
| MW | Molecular weight, including condensed species. |
| P | Pressure in bar. |
| Pr_eq | Equilibrium Prandtl number. |
| Pr_fr | Frozen Prandtl number. |
| T | Temperature in K. |
| conductivity_eq | Equilibrium thermal conductivity in \mu W/(cm·K). |
| conductivity_fr | Frozen thermal conductivity in \mu W/(cm·K). |
| converged | Convergence status of the solution. |
| cp | Equilibrium specific heat at constant pressure in kJ/(kg·K). |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). |
| cv | Equilibrium specific heat at constant volume in kJ/(kg·K). |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). |
| density | Density in kg/m³. |
| energy | Specific internal energy in kJ/kg. |
| enthalpy | Specific enthalpy in kJ/kg. |
| entropy | Specific entropy in kJ/(kg·K). |
| gamma_s | Isentropic specific heat ratio. |
| gibbs_energy | Specific Gibbs energy in kJ/kg. |
| last_error | int |
| ln_nj | Natural logarithm of gas phase species concentrations (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. |
| mass_fractions | Species mass fractions. |
| mole_fractions | Species mole fractions. |
| n | Total moles. |
| nj | Species concentrations (kg-mol of species j per kg of mixture). |
| viscosity | Viscosity in millipoise (mP). |
| volume | Specific volume in m³/kg. |

### EqDerivatives

Constructor: `EqDerivatives(solver, solution)`.

| Member | Signature / return convention |
|---|---|
| compute_derivatives | `compute_derivatives(self, check_closure_defect=False)` |
| compute_fd | `compute_fd(self, h, verbose=False, central=False)` |
| dG_dstate1 | See conventions above |
| dG_dstate1_fd | See conventions above |
| dG_dstate2 | See conventions above |
| dG_dstate2_fd | See conventions above |
| dG_dw0 | See conventions above |
| dG_dw0_fd | See conventions above |
| dH_dstate1 | See conventions above |
| dH_dstate1_fd | See conventions above |
| dH_dstate2 | See conventions above |
| dH_dstate2_fd | See conventions above |
| dH_dw0 | See conventions above |
| dH_dw0_fd | See conventions above |
| dS_dstate1 | See conventions above |
| dS_dstate1_fd | See conventions above |
| dS_dstate2 | See conventions above |
| dS_dstate2_fd | See conventions above |
| dS_dw0 | See conventions above |
| dS_dw0_fd | See conventions above |
| dT_dstate1 | See conventions above |
| dT_dstate1_fd | See conventions above |
| dT_dstate2 | See conventions above |
| dT_dstate2_fd | See conventions above |
| dT_dw0 | See conventions above |
| dT_dw0_fd | See conventions above |
| dU_dstate1 | See conventions above |
| dU_dstate1_fd | See conventions above |
| dU_dstate2 | See conventions above |
| dU_dstate2_fd | See conventions above |
| dU_dw0 | See conventions above |
| dU_dw0_fd | See conventions above |
| dn_dstate1 | See conventions above |
| dn_dstate1_fd | See conventions above |
| dn_dstate2 | See conventions above |
| dn_dstate2_fd | See conventions above |
| dn_dw0 | See conventions above |
| dn_dw0_fd | See conventions above |
| dnj_dstate1 | See conventions above |
| dnj_dstate1_fd | See conventions above |
| dnj_dstate2 | See conventions above |
| dnj_dstate2_fd | See conventions above |
| dnj_dw0 | See conventions above |
| dnj_dw0_fd | See conventions above |
| last_error | int |

### RocketSolution

Constructor: `RocketSolution(solver)`.

| Member | Signature / return convention |
|---|---|
| Isp | Specific impulse at each station in m/s. |
| Isp_vacuum | Vacuum specific impulse at each station in m/s. |
| M | Molecular weight (1/n) (gasses only) at each station. |
| MW | Molecular weight (including condensed species) at each station. |
| Mach | Mach number at each station. |
| P | Pressure at each station in bar. |
| Pr_eq | Equilibrium Prandtl number. |
| Pr_fr | Frozen Prandtl number. |
| T | Temperature at each station in K. |
| ae_at | Area ratio (Ae/At) at each station. |
| c_star | Characteristic velocity in m/s. |
| coefficient_of_thrust | Thrust coefficient at each station. |
| conductivity_eq | Equilibrium thermal conductivity in \mu W/(cm·K). |
| conductivity_fr | Frozen thermal conductivity in \mu W/(cm·K). |
| converged | Convergence status of the solution. |
| cp | Equilibrium specific heat at constant pressure in kJ/(kg·K). |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). |
| cv | Equilibrium specific heat at constant volume in kJ/(kg·K). |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). |
| density | Density at each station in kg/m³. |
| energy | Specific internal energy at each station in kJ/kg. |
| enthalpy | Specific enthalpy at each station in kJ/kg. |
| entropy | Specific entropy at each station in kJ/(kg·K). |
| gamma_s | Isentropic exponent \gamma_s at each station. |
| gibbs_energy | Specific Gibbs energy at each station in kJ/kg. |
| last_error | int |
| ln_nj | Natural logarithm of gas phase species concentrations at each station (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. |
| mass_fractions | Species mass fractions at each station. |
| mole_fractions | Species mole fractions at each station. |
| n | Total moles at each station. |
| nj | Species concentrations (kg-mol of species j per kg of mixture) at each station. |
| num_pts | Number of calculation points (stations). |
| sonic_velocity | Sonic velocity at each station in m/s. |
| viscosity | Viscosity in millipoise (mP). |
| volume | Specific volume at each station in m³/kg. |

### ShockSolution

Constructor: `ShockSolution(solver, reflected=True)`.

| Member | Signature / return convention |
|---|---|
| M | Molecular weight (1/n) at each shock state. |
| M21 | Molecular-weight ratio across the incident shock (M2/M1), not a Mach ratio. |
| M52 | Molecular-weight ratio across the reflected shock (M5/M2), not a Mach ratio. |
| MW | Molecular weight (including condensed species) at each shock state in kg/kmol. |
| Mach | Gas Mach number: ``velocity / sonic_velocity`` at each shock state. |
| P | Pressure at each shock state in bar. |
| P21 | Pressure ratio across the incident shock (P2/P1). |
| P52 | Pressure ratio across the reflected shock (P5/P2). |
| Pr_eq | Equilibrium Prandtl number at each shock state. |
| Pr_fr | Frozen Prandtl number at each shock state. |
| T | Temperature at each shock state in K. |
| T21 | Temperature ratio across the incident shock (T2/T1). |
| T52 | Temperature ratio across the reflected shock (T5/T2). |
| conductivity_eq | Equilibrium thermal conductivity at each shock state in \mu W/(cm·K). |
| conductivity_fr | Frozen thermal conductivity at each shock state in \mu W/(cm·K). |
| converged | Convergence status of the solution. |
| cp | Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K). |
| cp_eq | Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K). |
| cp_fr | Frozen specific heat at constant pressure at each shock state in kJ/(kg·K). |
| cv | Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K). |
| cv_eq | Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K). |
| cv_fr | Frozen specific heat at constant volume at each shock state in kJ/(kg·K). |
| density | Density at each shock state in kg/m³. |
| energy | Specific internal energy at each shock state in kJ/kg. |
| enthalpy | Specific enthalpy at each shock state in kJ/kg. |
| entropy | Specific entropy at each shock state in kJ/(kg·K). |
| gamma_s | Isentropic exponent \gamma_s at each shock state. |
| gibbs_energy | Gibbs energy at each shock state in kJ/kg. |
| last_error | int |
| ln_nj | Natural logarithm of gas phase species concentrations at each station (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. |
| mass_fractions | Species mass fractions at each station. |
| mole_fractions | Species mole fractions at each station. |
| n | Total moles at each station. |
| nj | Species concentrations (kg-mol of species j per kg of mixture) at each index. |
| rho12 | Density ratio across the incident shock (ρ1/ρ2). |
| rho52 | Density ratio across the reflected shock (ρ5/ρ2). |
| sonic_velocity | Reported sound speed in m/s: ``sqrt(R * T * gamma_s / M)``. |
| u5_p_v2 | Upstream gas speed relative to the reflected shock in m/s. |
| v2 | Station-2 gas speed in the wall/unshocked-gas frame in m/s. |
| velocity | Gas speed in m/s, in each station's shock frame. |
| viscosity | Viscosity at each shock state in millipoise (mP). |
| volume | Specific volume at each shock state in m³/kg. |

### DetonationSolution

Constructor: `DetonationSolution(solver)`.

| Member | Signature / return convention |
|---|---|
| H1 | Initial specific enthalpy in kJ/kg. |
| M | Molecular weight (1/n) (gasses only) of detonation products in kg/kmol. |
| M1 | Initial molecular weight (1/n). |
| MW | Molecular weight (including condensed species) of detonation products in kg/kmol. |
| M_M1 | Molecular weight ratio M/M1. |
| Mach | Detonation Mach number. |
| P | Detonation pressure in bar. |
| P1 | Initial pressure in bar. |
| P_P1 | Pressure ratio P/P1. |
| Pr_eq | Equilibrium Prandtl number. |
| Pr_fr | Frozen Prandtl number. |
| T | Detonation temperature in K. |
| T1 | Initial temperature in K. |
| T_T1 | Temperature ratio T/T1. |
| conductivity_eq | Equilibrium thermal conductivity in \mu W/(cm·K). |
| conductivity_fr | Frozen thermal conductivity in \mu W/(cm·K). |
| converged | Convergence status of the solution. |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). |
| density | Detonation density in kg/m³. |
| energy | Detonation specific internal energy in kJ/kg. |
| enthalpy | Detonation specific enthalpy in kJ/kg. |
| entropy | Detonation specific entropy in kJ/(kg·K). |
| gamma1 | Initial isentropic exponent \gamma_1. |
| gamma_s | Detonation isentropic exponent \gamma_s. |
| gibbs_energy | Detonation Gibbs energy in kJ/kg. |
| last_error | int |
| ln_nj | Natural logarithm of gas phase species concentrations (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. |
| mass_fractions | Species mass fractions. |
| mole_fractions | Species mole fractions. |
| n | Total moles. |
| nj | Species concentrations (kg-mol of species j per kg of mixture). |
| rho_rho1 | Density ratio ρ/ρ1. |
| sonic_velocity | Detonation sonic velocity in m/s. |
| sonic_velocity1 | Initial sonic velocity in m/s. |
| velocity | Detonation velocity in m/s. |
| viscosity | Viscosity in millipoise (mP). |

### RocketSolver

Constructor: `RocketSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`.

| Member | Signature / return convention |
|---|---|
| num_condensed | Number of condensed phase species. |
| num_elements | Number of chemical elements. |
| num_gas | Number of gas phase species. |
| num_products | Number of product species. |
| num_reactants | Number of reactant species. |
| solve | `solve(self, soln, weights, pc, pi_p=None, subar=None, supar=None, iac=True, n_frz=None, hc=None, tc=None, mdot=None, ac_at=None, tc_est=None)` |

### ShockSolver

Constructor: `ShockSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`.

| Member | Signature / return convention |
|---|---|
| num_condensed | Number of condensed phase species. |
| num_elements | Number of chemical elements. |
| num_gas | Number of gas phase species. |
| num_products | Number of product species. |
| num_reactants | Number of reactant species. |
| solve | `solve(self, soln, weights, T0, p0, u1=None, Mach1=None, reflected=True, incident_frozen=False, reflected_frozen=False)` |

### DetonationSolver

Constructor: `DetonationSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`.

| Member | Signature / return convention |
|---|---|
| num_condensed | Number of condensed phase species. |
| num_elements | Number of chemical elements. |
| num_gas | Number of gas phase species. |
| num_products | Number of product species. |
| num_reactants | Number of reactant species. |
| solve | `solve(self, soln, weights, T1, p1, frozen=False)` |

## Root functions

| Member | Signature |
|---|---|
| init | `init(path=None, trans=True)` |
| init_thermo | `init_thermo(thermofile)` |
| init_trans | `init_trans(transfile)` |
| is_initialized | `is_initialized()` |
| set_log_level | `set_log_level(level)` |
| eq_solve | `eq_solve(eq_type, reactants, T: 'Optional[double]' = None, H: 'Optional[double]' = None, S: 'Optional[double]' = None, U: 'Optional[double]' = None, P: 'Optional[double]' = None, V: 'Optional[double]' = None, T_reac: 'Optional[list | np.ndarray | double]' = None, fuel_amounts: 'Optional[list | np.ndarray]' = None, oxid_amounts: 'Optional[list | np.ndarray]' = None, moles=False, of_ratio: 'Optional[list | np.ndarray | double]' = None, phi: 'Optional[list | np.ndarray | double]' = None, r_eq: 'Optional[list | np.ndarray | double]' = None, pct_fuel: 'Optional[list | np.ndarray | double]' = None, only: 'Optional[list]' = None, omit: 'Optional[list]' = None, insert: 'Optional[list]' = None, trace: 'Optional[double]' = None, transport=False, ions=False)` |
| lib_version | `lib_version()` |
| lib_version_major | `lib_version_major()` |
| lib_version_minor | `lib_version_minor()` |
| lib_version_patch | `lib_version_patch()` |

__version__ is the Python distribution version; lib_version* identify the actual
backend. The adapter reports its own Rust version, never claims to be NASA 3.3.4.

## Constants

All integer constants are dimensionless selectors or status codes.

| Member | Value |
|---|---|
| SUCCESS | 0 |
| INVALID_FILENAME | 1 |
| INVALID_PROPERTY_TYPE | 2 |
| INVALID_EQUILIBRIUM_TYPE | 3 |
| NOT_CONVERGED | 8 |
| LAST_VALID_SOLUTION | 9 |
| FORTRAN_ABORT | 10 |
| LOG_CRITICAL | 50 |
| LOG_ERROR | 40 |
| LOG_WARNING | 30 |
| LOG_INFO | 20 |
| LOG_DEBUG | 10 |
| LOG_NONE | 0 |
| TP | 0 |
| HP | 1 |
| SP | 2 |
| TV | 3 |
| UV | 4 |
| SV | 5 |
| NUM_REACTANTS | 0 |
| NUM_PRODUCTS | 1 |
| NUM_GAS | 2 |
| NUM_CONDENSED | 3 |
| NUM_ELEMENTS | 4 |
| MAX_EQUATIONS | 5 |
| TEMPERATURE | 0 |
| PRESSURE | 1 |
| VOLUME | 2 |
| DENSITY | 3 |
| ENTHALPY | 6 |
| ENERGY | 7 |
| ENTROPY | 8 |
| GIBBS_ENERGY | 9 |
| GAMMA_S | 10 |
| FROZEN_CP | 11 |
| FROZEN_CV | 12 |
| EQUILIBRIUM_CP | 13 |
| EQUILIBRIUM_CV | 14 |
| VISCOSITY | 15 |
| FROZEN_CONDUCTIVITY | 16 |
| EQUILIBRIUM_CONDUCTIVITY | 17 |
| FROZEN_PRANDTL | 18 |
| EQUILIBRIUM_PRANDTL | 19 |
| ROCKET_TEMPERATURE | 0 |
| ROCKET_PRESSURE | 1 |
| ROCKET_VOLUME | 2 |
| ROCKET_DENSITY | 3 |
| ROCKET_M | 4 |
| ROCKET_MW | 5 |
| ROCKET_ENTHALPY | 6 |
| ROCKET_ENERGY | 7 |
| ROCKET_ENTROPY | 8 |
| ROCKET_GIBBS_ENERGY | 9 |
| ROCKET_GAMMA_S | 10 |
| ROCKET_FROZEN_CP | 11 |
| ROCKET_FROZEN_CV | 12 |
| ROCKET_EQUILIBRIUM_CP | 13 |
| ROCKET_EQUILIBRIUM_CV | 14 |
| MACH | 15 |
| SONIC_VELOCITY | 16 |
| AE_AT | 17 |
| C_STAR | 18 |
| COEFFICIENT_OF_THRUST | 19 |
| ISP | 20 |
| ISP_VACUUM | 21 |
| ROCKET_VISCOSITY | 22 |
| ROCKET_FROZEN_CONDUCTIVITY | 23 |
| ROCKET_EQUILIBRIUM_CONDUCTIVITY | 24 |
| ROCKET_FROZEN_PRANDTL | 25 |
| ROCKET_EQUILIBRIUM_PRANDTL | 26 |
| SHOCK_TEMPERATURE | 0 |
| SHOCK_PRESSURE | 1 |
| SHOCK_VELOCITY | 2 |
| SHOCK_MACH | 3 |
| SHOCK_SONIC_VELOCITY | 4 |
| SHOCK_RHO12 | 5 |
| SHOCK_RHO52 | 6 |
| SHOCK_P21 | 7 |
| SHOCK_P52 | 8 |
| SHOCK_T21 | 9 |
| SHOCK_T52 | 10 |
| SHOCK_M21 | 11 |
| SHOCK_M52 | 12 |
| SHOCK_V2 | 13 |
| SHOCK_U5_P_V2 | 14 |
| SHOCK_VOLUME | 15 |
| SHOCK_DENSITY | 16 |
| SHOCK_M | 17 |
| SHOCK_MW | 18 |
| SHOCK_ENTHALPY | 19 |
| SHOCK_ENERGY | 20 |
| SHOCK_ENTROPY | 21 |
| SHOCK_GIBBS_ENERGY | 22 |
| SHOCK_GAMMA_S | 23 |
| SHOCK_FROZEN_CP | 24 |
| SHOCK_FROZEN_CV | 25 |
| SHOCK_EQUILIBRIUM_CP | 26 |
| SHOCK_EQUILIBRIUM_CV | 27 |
| SHOCK_VISCOSITY | 28 |
| SHOCK_FROZEN_CONDUCTIVITY | 29 |
| SHOCK_EQUILIBRIUM_CONDUCTIVITY | 30 |
| SHOCK_FROZEN_PRANDTL | 31 |
| SHOCK_EQUILIBRIUM_PRANDTL | 32 |
| DETONATION_P1 | 0 |
| DETONATION_T1 | 1 |
| DETONATION_H1 | 2 |
| DETONATION_M1 | 3 |
| DETONATION_GAMMA1 | 4 |
| DETONATION_V_SONIC1 | 5 |
| DETONATION_PRESSURE | 6 |
| DETONATION_TEMPERATURE | 7 |
| DETONATION_DENSITY | 8 |
| DETONATION_ENTHALPY | 9 |
| DETONATION_ENERGY | 10 |
| DETONATION_GIBBS_ENERGY | 11 |
| DETONATION_ENTROPY | 12 |
| DETONATION_MACH | 13 |
| DETONATION_VELOCITY | 14 |
| DETONATION_SONIC_VELOCITY | 15 |
| DETONATION_GAMMA | 16 |
| DETONATION_P_P1 | 17 |
| DETONATION_T_T1 | 18 |
| DETONATION_M_M1 | 19 |
| DETONATION_RHO_RHO1 | 20 |
| DETONATION_FROZEN_CP | 21 |
| DETONATION_FROZEN_CV | 22 |
| DETONATION_EQUILIBRIUM_CP | 23 |
| DETONATION_EQUILIBRIUM_CV | 24 |
| DETONATION_M | 25 |
| DETONATION_MW | 26 |
| DETONATION_VISCOSITY | 27 |
| DETONATION_FROZEN_CONDUCTIVITY | 28 |
| DETONATION_EQUILIBRIUM_CONDUCTIVITY | 29 |
| DETONATION_FROZEN_PRANDTL | 30 |
| DETONATION_EQUILIBRIUM_PRANDTL | 31 |
| R / constants.R | 8314.51 J/(kmol K) |

## Unit helpers

Each takes one argument, value; no default; returns the same scalar/NumPy shape.

| Member | Conversion |
|---|---|
| celsius_to_kelvin(value) | Convert °C to K. |
| kelvin_to_celsius(value) | Convert K to °C. |
| fahrenheit_to_kelvin(value) | Convert °F to K. |
| kelvin_to_fahrenheit(value) | Convert K to °F. |
| rankine_to_kelvin(value) | Convert °R to K. |
| kelvin_to_rankine(value) | Convert K to °R. |
| atm_to_bar(value) | Convert atmospheres to bar. |
| bar_to_atm(value) | Convert bar to atmospheres. |
| psi_to_bar(value) | Convert psi to bar. |
| bar_to_psi(value) | Convert bar to psi. |
| mmhg_to_bar(value) | Convert mmHg to bar. |
| bar_to_mmhg(value) | Convert bar to mmHg. |
| g_per_cm3_to_kg_per_m3(value) | Convert g/cm³ to kg/m³. |
| kg_per_m3_to_g_per_cm3(value) | Convert kg/m³ to g/cm³. |
| cm3_per_g_to_m3_per_kg(value) | Convert cm³/g to m³/kg. |
| m3_per_kg_to_cm3_per_g(value) | Convert m³/kg to cm³/g. |
| m_per_s_to_ft_per_s(value) | Convert m/s to ft/s. |
| ft_per_s_to_m_per_s(value) | Convert ft/s to m/s. |
| joule_to_cal(value) | Convert joules to calories. |
| cal_to_joule(value) | Convert calories to joules. |
| kj_per_kg_to_cal_per_g(value) | Convert kJ/kg to cal/g. |
| cal_per_g_to_kj_per_kg(value) | Convert cal/g to kJ/kg. |

## Other namespaces

cea.lib re-exports the same library API. The optional MATLAB module exposes
eq_solve, rocket_solve, shock_solve and detonation_solve, returning flattened
SimpleNamespace objects; its Matlab*Solution stub types carry the same fields
as their corresponding solutions. MATLAB wrappers, sample plotting and file
outputs are outside this adapter. Access raises a named NotImplementedError.
The native __all__ also accidentally exports os and importlib_resources;
these helper modules are not re-exported.
