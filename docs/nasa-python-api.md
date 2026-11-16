# NASA Python 3.3.4 API contract

Inspected the PyPI wheel and Python/Cython sources at NASA commit
[2f79a647](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da/source/bind/python)
on 2026-11-17. This inventory records public API facts, not copied implementation.
NASA is Apache-2.0; see [provenance](provenance.md) and the root NOTICE.

Install `cea-keplr-rs[nasa]` to use this adapter. NumPy is optional for the
base package: `import cea_rs` and its SI API work without it. Importing
`cea_rs.nasa` without NumPy raises ImportError naming that extra.

The largest compared relative difference, including named exceptions, is
**1.73451869e-4**: RP-1311 example 12's CO mass fraction at area 75, within
its 1.8e-4 bound. The worst other valid compared output is 7.0467525e-5.
Fixed-constraint finite-difference noise and invalid native states have the
separate bounds described below.

## Conventions

- Solver pressure inputs and solution pressures: bar. Temperatures: K.
- EqSolver state1: T for TP/TV, H/R for HP, S/R for SP/SV, U/R for UV,
  where H/U are J/kg, S is J/(kg K), and R=8314.51 J/(kmol K).
  State2 is bar for TP/HP/SP and m³/kg for TV/UV/SV.
- Rocket hc uses H/R despite its docstring saying (kJ/kg)/R. mdot is kg/(s m²).
- Solutions: H/U/G in kJ/kg, Cp/Cv/S in kJ/(kg K), volume m³/kg, density kg/m³,
  M/MW kg/kmol, velocities/c-star/Isp in m/s, viscosity millipoise,
  conductivity mW/(cm K) (native docstrings incorrectly say µW); ratios, fractions and Prandtl numbers are dimensionless.
- Mixture.calc_property: H/U/G J/kg, scalar-temperature S J/(kg K), Cp/Cv J/(kg K).
  In 3.3.4 array-temperature S is divided by 1000; array-temperature G retains
  the supplied weight scale and divides bar pressure by 1e5; V/DENSITY retain legacy ×1000/×0.001 relative to SI.
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

Constructor: `Reactant(name, formula=None, molecular_weight=None, enthalpy=None, enthalpy_units=None, temperature=None)`. Default constructor implemented. Overrides of fitted named reactants raise NotImplementedError; custom reference reactants and reference-only overrides are supported.

NASA does not expose a public Python phase property. Internally, the adapter
uses compiled species phases and retains parsed reference-reactant phases.
For example, NASA 3.3.4's bundled thermo.lib has Co(b) interval phases 3–6;
its ThermoDB loader retains the final value, 6, matching thermo_species.
The adapter no longer scans thermo.inp on each phase query.

| Member | Signature / return convention | Status |
|---|---|---|
| enthalpy | See conventions above | Implemented |
| enthalpy_units | See conventions above | Implemented |
| formula | See conventions above | Implemented |
| get_valid_temperature_range | `get_valid_temperature_range(self)` | Implemented |
| molecular_weight | See conventions above | Implemented |
| name | See conventions above | Implemented |
| temperature | See conventions above | Implemented |

### Mixture

Constructor: `Mixture(species=None, products_from_reactants=False, omit=[], ions=False)`. Default constructor implemented.

NASA accepts `omit=None` as an empty omission list. The adapter keeps NASA's
visible `omit=[]` default and never mutates that list.

| Member | Signature / return convention | Status |
|---|---|---|
| calc_property | `calc_property(self, prop_type, weights, temperature: 'float \| list \| np.ndarray', pressure: 'Optional[float]' = None)` | Implemented for available fits; reference-only Cp/S raises NotImplementedError |
| chem_eq_ratio_to_of_ratio | `chem_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, chem_eq_ratio)` | Implemented |
| moles_to_weights | `moles_to_weights(self, moles)` | Implemented |
| num_species | int | Implemented |
| of_ratio_to_weights | `of_ratio_to_weights(self, oxidant_weights, fuel_weights, of_ratio)` | Implemented |
| species_names | Get list of species names in the mixture. | Implemented |
| weight_eq_ratio_to_of_ratio | `weight_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, weight_eq_ratio)` | Implemented |
| weights_to_moles | `weights_to_moles(self, weights)` | Implemented |

### EqSolver

Constructor: `EqSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`. Default constructor implemented. smooth_truncation=True or truncation_width>0 raises NotImplementedError.

| Member | Signature / return convention | Status |
|---|---|---|
| max_equations | Maximum number of equations in system; num_elements + num_condensed + (1 for TP/TV problems, 2 otherwise). | Implemented |
| num_condensed | Number of condensed phase species. | Implemented |
| num_elements | Number of chemical elements. | Implemented |
| num_gas | Number of gas phase species. | Implemented |
| num_products | Number of product species. | Implemented |
| num_reactants | Number of reactant species. | Implemented |
| solve | `solve(self, soln, eq_type, state1, state2, amounts)` | Implemented |

### EqSolution

Constructor: `EqSolution(solver, T_init=None, nj_init=None)`. Default constructor implemented. T_init or nj_init raises NotImplementedError.

| Member | Signature / return convention | Status |
|---|---|---|
| M | Molecular weight, gasses only (1/n). | Implemented |
| MW | Molecular weight, including condensed species. | Implemented |
| P | Pressure in bar. | Implemented |
| Pr_eq | Equilibrium Prandtl number. | Implemented |
| Pr_fr | Frozen Prandtl number. | Implemented |
| T | Temperature in K. | Implemented |
| conductivity_eq | Equilibrium thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| conductivity_fr | Frozen thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| converged | Convergence status of the solution. | Implemented; failure clears old results and raises instead of retaining invalid native values |
| cp | Equilibrium specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cv | Equilibrium specific heat at constant volume in kJ/(kg·K). | Implemented |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). | Implemented |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). | Implemented |
| density | Density in kg/m³. | Implemented |
| energy | Specific internal energy in kJ/kg. | Implemented |
| enthalpy | Specific enthalpy in kJ/kg. | Implemented |
| entropy | Specific entropy in kJ/(kg·K). | Implemented |
| gamma_s | Isentropic specific heat ratio. | Implemented |
| gibbs_energy | Specific Gibbs energy in kJ/kg. | Implemented |
| last_error | int | Implemented |
| ln_nj | Natural logarithm of gas phase species concentrations (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. | Implemented |
| mass_fractions | Species mass fractions. | Implemented |
| mole_fractions | Species mole fractions. | Implemented |
| n | Total moles. | Implemented |
| nj | Species concentrations (kg-mol of species j per kg of mixture). | Implemented |
| viscosity | Viscosity in millipoise (mP). | Implemented |
| volume | Specific volume in m³/kg. | Implemented |

### EqDerivatives

Constructor: `EqDerivatives(solver, solution)`. Default constructor implemented.

NASA's `check_closure_defect=True` evaluates the internal residual
`J * du/dx + Rx`; it does not change the returned derivatives. Its Python
docstring mentions logging, but the Fortran check stores the residual without
printing unless verbose output is requested separately. Rust already evaluates
the same residual and rejects a scaled defect above 1e-10 on every analytic
solve. The adapter's explicit flag verifies that the core supplied a finite
defect within that bound, raising ConvergenceError otherwise. The default
retains the core's existing validation.

| Member | Signature / return convention | Status |
|---|---|---|
| compute_derivatives | `compute_derivatives(self, check_closure_defect=False)` | Implemented |
| compute_fd | `compute_fd(self, h, verbose=False, central=False)` | Implemented |
| dG_dstate1 | See conventions above | Implemented |
| dG_dstate1_fd | See conventions above | Implemented |
| dG_dstate2 | See conventions above | Implemented |
| dG_dstate2_fd | See conventions above | Implemented |
| dG_dw0 | See conventions above | Implemented |
| dG_dw0_fd | See conventions above | Implemented |
| dH_dstate1 | See conventions above | Implemented |
| dH_dstate1_fd | See conventions above | Implemented |
| dH_dstate2 | See conventions above | Implemented |
| dH_dstate2_fd | See conventions above | Implemented |
| dH_dw0 | See conventions above | Implemented |
| dH_dw0_fd | See conventions above | Implemented |
| dS_dstate1 | See conventions above | Implemented |
| dS_dstate1_fd | See conventions above | Implemented |
| dS_dstate2 | See conventions above | Implemented |
| dS_dstate2_fd | See conventions above | Implemented |
| dS_dw0 | See conventions above | Implemented |
| dS_dw0_fd | See conventions above | Implemented |
| dT_dstate1 | See conventions above | Implemented |
| dT_dstate1_fd | See conventions above | Implemented |
| dT_dstate2 | See conventions above | Implemented |
| dT_dstate2_fd | See conventions above | Implemented |
| dT_dw0 | See conventions above | Implemented |
| dT_dw0_fd | See conventions above | Implemented |
| dU_dstate1 | See conventions above | Implemented |
| dU_dstate1_fd | See conventions above | Implemented |
| dU_dstate2 | See conventions above | Implemented |
| dU_dstate2_fd | See conventions above | Implemented |
| dU_dw0 | See conventions above | Implemented |
| dU_dw0_fd | See conventions above | Implemented |
| dn_dstate1 | See conventions above | Implemented |
| dn_dstate1_fd | See conventions above | Implemented |
| dn_dstate2 | See conventions above | Implemented |
| dn_dstate2_fd | See conventions above | Implemented |
| dn_dw0 | See conventions above | Implemented |
| dn_dw0_fd | See conventions above | Implemented |
| dnj_dstate1 | See conventions above | Implemented |
| dnj_dstate1_fd | See conventions above | Implemented |
| dnj_dstate2 | See conventions above | Implemented |
| dnj_dstate2_fd | See conventions above | Implemented |
| dnj_dw0 | See conventions above | Implemented |
| dnj_dw0_fd | See conventions above | Implemented |
| last_error | int | Implemented |

### RocketSolution

Constructor: `RocketSolution(solver)`. Default constructor implemented.

| Member | Signature / return convention | Status |
|---|---|---|
| Isp | Specific impulse at each station in m/s. | Implemented |
| Isp_vacuum | Vacuum specific impulse at each station in m/s. | Implemented |
| M | Molecular weight (1/n) (gasses only) at each station. | Implemented |
| MW | Molecular weight (including condensed species) at each station. | Implemented |
| Mach | Mach number at each station. | Implemented |
| P | Pressure at each station in bar. | Implemented |
| Pr_eq | Equilibrium Prandtl number. | Implemented |
| Pr_fr | Frozen Prandtl number. | Implemented |
| T | Temperature at each station in K. | Implemented |
| ae_at | Area ratio (Ae/At) at each station. | Implemented |
| c_star | Characteristic velocity in m/s. | Implemented |
| coefficient_of_thrust | Thrust coefficient at each station. | Implemented |
| conductivity_eq | Equilibrium thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| conductivity_fr | Frozen thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| converged | Convergence status of the solution. | Implemented; failure clears old results and raises instead of retaining invalid native values |
| cp | Equilibrium specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cv | Equilibrium specific heat at constant volume in kJ/(kg·K). | Differs after freezing: current Cv instead of retained/zero native fields |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). | Differs after freezing: current Cv instead of retained/zero native fields |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). | Differs after freezing: current Cv instead of retained/zero native fields |
| density | Density at each station in kg/m³. | Implemented |
| energy | Specific internal energy at each station in kJ/kg. | Implemented |
| enthalpy | Specific enthalpy at each station in kJ/kg. | Implemented |
| entropy | Specific entropy at each station in kJ/(kg·K). | Implemented |
| gamma_s | Isentropic exponent \gamma_s at each station. | Implemented |
| gibbs_energy | Specific Gibbs energy at each station in kJ/kg. | Implemented |
| last_error | int | Implemented |
| ln_nj | Natural logarithm of gas phase species concentrations at each station (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. | NotImplementedError when any station lacks retained equilibrium logs (all shocks include station 1) |
| mass_fractions | Species mass fractions at each station. | Implemented |
| mole_fractions | Species mole fractions at each station. | Implemented |
| n | Total moles at each station. | Implemented |
| nj | Species concentrations (kg-mol of species j per kg of mixture) at each station. | Implemented |
| num_pts | Number of calculation points (stations). | Implemented |
| sonic_velocity | Sonic velocity at each station in m/s. | Implemented |
| viscosity | Viscosity in millipoise (mP). | Implemented |
| volume | Specific volume at each station in m³/kg. | Implemented |

### ShockSolution

Constructor: `ShockSolution(solver, reflected=True)`. Default constructor implemented.

| Member | Signature / return convention | Status |
|---|---|---|
| M | Molecular weight (1/n) at each shock state. | Implemented |
| M21 | Molecular-weight ratio across the incident shock (M2/M1), not a Mach ratio. | Implemented |
| M52 | Molecular-weight ratio across the reflected shock (M5/M2), not a Mach ratio. | Implemented |
| MW | Molecular weight (including condensed species) at each shock state in kg/kmol. | Implemented |
| Mach | Gas Mach number: ``velocity / sonic_velocity`` at each shock state. | Implemented |
| P | Pressure at each shock state in bar. | Implemented |
| P21 | Pressure ratio across the incident shock (P2/P1). | Implemented |
| P52 | Pressure ratio across the reflected shock (P5/P2). | Implemented |
| Pr_eq | Equilibrium Prandtl number at each shock state. | Implemented |
| Pr_fr | Frozen Prandtl number at each shock state. | Implemented |
| T | Temperature at each shock state in K. | Implemented |
| T21 | Temperature ratio across the incident shock (T2/T1). | Implemented |
| T52 | Temperature ratio across the reflected shock (T5/T2). | Implemented |
| conductivity_eq | Equilibrium thermal conductivity at each shock state in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| conductivity_fr | Frozen thermal conductivity at each shock state in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| converged | Convergence status of the solution. | Implemented; failure clears old results and raises instead of retaining invalid native values |
| cp | Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| cp_eq | Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| cp_fr | Frozen specific heat at constant pressure at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| cv | Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| cv_eq | Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| cv_fr | Frozen specific heat at constant volume at each shock state in kJ/(kg·K). | Differs at frozen station 5: computes defined Cp/Cv instead of uninitialized native fields |
| density | Density at each shock state in kg/m³. | Implemented |
| energy | Specific internal energy at each shock state in kJ/kg. | Implemented |
| enthalpy | Specific enthalpy at each shock state in kJ/kg. | Implemented |
| entropy | Specific entropy at each shock state in kJ/(kg·K). | Implemented |
| gamma_s | Isentropic exponent \gamma_s at each shock state. | Implemented |
| gibbs_energy | Gibbs energy at each shock state in kJ/kg. | Implemented |
| last_error | int | Implemented |
| ln_nj | Natural logarithm of gas phase species concentrations at each station (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. | NotImplementedError when any station lacks retained equilibrium logs (all shocks include station 1) |
| mass_fractions | Species mass fractions at each station. | Implemented |
| mole_fractions | Species mole fractions at each station. | Implemented |
| n | Total moles at each station. | Differs: returns documented arrays; native getter raises AttributeError |
| nj | Species concentrations (kg-mol of species j per kg of mixture) at each index. | Differs: returns documented arrays; native getter raises AttributeError |
| rho12 | Density ratio across the incident shock (ρ1/ρ2). | Implemented |
| rho52 | Density ratio across the reflected shock (ρ5/ρ2). | Implemented |
| sonic_velocity | Reported sound speed in m/s: ``sqrt(R * T * gamma_s / M)``. | Implemented |
| u5_p_v2 | Upstream gas speed relative to the reflected shock in m/s. | Implemented |
| v2 | Station-2 gas speed in the wall/unshocked-gas frame in m/s. | Implemented |
| velocity | Gas speed in m/s, in each station's shock frame. | Implemented |
| viscosity | Viscosity at each shock state in millipoise (mP). | Implemented |
| volume | Specific volume at each shock state in m³/kg. | Implemented |

### DetonationSolution

Constructor: `DetonationSolution(solver)`. Default constructor implemented.

| Member | Signature / return convention | Status |
|---|---|---|
| H1 | Initial specific enthalpy in kJ/kg. | Implemented |
| M | Molecular weight (1/n) (gasses only) of detonation products in kg/kmol. | Implemented |
| M1 | Initial molecular weight (1/n). | Implemented |
| MW | Molecular weight (including condensed species) of detonation products in kg/kmol. | Implemented |
| M_M1 | Molecular weight ratio M/M1. | Implemented |
| Mach | Detonation Mach number. | Implemented |
| P | Detonation pressure in bar. | Implemented |
| P1 | Initial pressure in bar. | Implemented |
| P_P1 | Pressure ratio P/P1. | Implemented |
| Pr_eq | Equilibrium Prandtl number. | Implemented |
| Pr_fr | Frozen Prandtl number. | Implemented |
| T | Detonation temperature in K. | Implemented |
| T1 | Initial temperature in K. | Implemented |
| T_T1 | Temperature ratio T/T1. | Implemented |
| conductivity_eq | Equilibrium thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| conductivity_fr | Frozen thermal conductivity in mW/(cm·K) (native docstring incorrectly says µW). | Implemented |
| converged | Convergence status of the solution. | Implemented; failure clears old results and raises instead of retaining invalid native values |
| cp_eq | Equilibrium specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cp_fr | Frozen specific heat at constant pressure in kJ/(kg·K). | Implemented |
| cv_eq | Equilibrium specific heat at constant volume in kJ/(kg·K). | Implemented |
| cv_fr | Frozen specific heat at constant volume in kJ/(kg·K). | Implemented |
| density | Detonation density in kg/m³. | Implemented |
| energy | Detonation specific internal energy in kJ/kg. | Implemented |
| enthalpy | Detonation specific enthalpy in kJ/kg. | Implemented |
| entropy | Detonation specific entropy in kJ/(kg·K). | Implemented |
| gamma1 | Initial isentropic exponent \gamma_1. | Implemented |
| gamma_s | Detonation isentropic exponent \gamma_s. | Implemented |
| gibbs_energy | Detonation Gibbs energy in kJ/kg. | Implemented |
| last_error | int | Implemented |
| ln_nj | Natural logarithm of gas phase species concentrations (gasses only). NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts. | Implemented |
| mass_fractions | Species mass fractions. | Implemented |
| mole_fractions | Species mole fractions. | Implemented |
| n | Total moles. | Implemented |
| nj | Species concentrations (kg-mol of species j per kg of mixture). | Implemented |
| rho_rho1 | Density ratio ρ/ρ1. | Implemented |
| sonic_velocity | Detonation sonic velocity in m/s. | Implemented |
| sonic_velocity1 | Initial sonic velocity in m/s. | Implemented |
| velocity | Detonation velocity in m/s. | Implemented |
| viscosity | Viscosity in millipoise (mP). | Implemented |

### RocketSolver

Constructor: `RocketSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`. Default constructor implemented. smooth_truncation=True or truncation_width>0 raises NotImplementedError.

| Member | Signature / return convention | Status |
|---|---|---|
| num_condensed | Number of condensed phase species. | Implemented |
| num_elements | Number of chemical elements. | Implemented |
| num_gas | Number of gas phase species. | Implemented |
| num_products | Number of product species. | Implemented |
| num_reactants | Number of reactant species. | Implemented |
| solve | `solve(self, soln, weights, pc, pi_p=None, subar=None, supar=None, iac=True, n_frz=None, hc=None, tc=None, mdot=None, ac_at=None, tc_est=None)` | Implemented |

### ShockSolver

Constructor: `ShockSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`. Default constructor implemented. smooth_truncation=True or truncation_width>0 raises NotImplementedError.

| Member | Signature / return convention | Status |
|---|---|---|
| num_condensed | Number of condensed phase species. | Implemented |
| num_elements | Number of chemical elements. | Implemented |
| num_gas | Number of gas phase species. | Implemented |
| num_products | Number of product species. | Implemented |
| num_reactants | Number of reactant species. | Implemented |
| solve | `solve(self, soln, weights, T0, p0, u1=None, Mach1=None, reflected=True, incident_frozen=False, reflected_frozen=False)` | Implemented |

### DetonationSolver

Constructor: `DetonationSolver(products, *, reactants=<products>, transport=False, ions=False, trace=-1.0, insert=[], smooth_truncation=False, truncation_width=-1.0)`. Default constructor implemented. smooth_truncation=True or truncation_width>0 raises NotImplementedError.

| Member | Signature / return convention | Status |
|---|---|---|
| num_condensed | Number of condensed phase species. | Implemented |
| num_elements | Number of chemical elements. | Implemented |
| num_gas | Number of gas phase species. | Implemented |
| num_products | Number of product species. | Implemented |
| num_reactants | Number of reactant species. | Implemented |
| solve | `solve(self, soln, weights, T1, p1, frozen=False)` | Implemented |

## Root functions

| Member | Signature | Status |
|---|---|---|
| init | `init(path=None, trans=True)` | Implemented for path=None; custom path raises NotImplementedError; trans is stateless |
| init_thermo | `init_thermo(thermofile)` | NotImplementedError: embedded databases only |
| init_trans | `init_trans(transfile)` | NotImplementedError: embedded databases only |
| is_initialized | `is_initialized()` | Implemented |
| set_log_level | `set_log_level(level)` | Silent levels accepted; debug/info raises NotImplementedError (no native logger) |
| eq_solve | `eq_solve(eq_type, reactants, T: 'Optional[double]' = None, H: 'Optional[double]' = None, S: 'Optional[double]' = None, U: 'Optional[double]' = None, P: 'Optional[double]' = None, V: 'Optional[double]' = None, T_reac: 'Optional[list \| np.ndarray \| double]' = None, fuel_amounts: 'Optional[list \| np.ndarray]' = None, oxid_amounts: 'Optional[list \| np.ndarray]' = None, moles=False, of_ratio: 'Optional[list \| np.ndarray \| double]' = None, phi: 'Optional[list \| np.ndarray \| double]' = None, r_eq: 'Optional[list \| np.ndarray \| double]' = None, pct_fuel: 'Optional[list \| np.ndarray \| double]' = None, only: 'Optional[list]' = None, omit: 'Optional[list]' = None, insert: 'Optional[list]' = None, trace: 'Optional[double]' = None, transport=False, ions=False)` | Implemented legacy root helper; preserves actual H/S/U-to-state1 behavior, despite J-based docstrings |
| lib_version | `lib_version()` | Differs: reports actual Rust version |
| lib_version_major | `lib_version_major()` | Differs: reports actual Rust version |
| lib_version_minor | `lib_version_minor()` | Differs: reports actual Rust version |
| lib_version_patch | `lib_version_patch()` | Differs: reports actual Rust version |

__version__ is the Python distribution version; lib_version* identify the actual
backend. The adapter reports its own Rust version, never claims to be NASA 3.3.4.

## Constants

All integer constants are dimensionless selectors or status codes.

| Member | Value | Status |
|---|---|---|
| SUCCESS | 0 | Implemented |
| INVALID_FILENAME | 1 | Implemented |
| INVALID_PROPERTY_TYPE | 2 | Implemented |
| INVALID_EQUILIBRIUM_TYPE | 3 | Implemented |
| NOT_CONVERGED | 8 | Implemented |
| LAST_VALID_SOLUTION | 9 | Implemented |
| FORTRAN_ABORT | 10 | Implemented |
| LOG_CRITICAL | 50 | Implemented |
| LOG_ERROR | 40 | Implemented |
| LOG_WARNING | 30 | Implemented |
| LOG_INFO | 20 | Implemented |
| LOG_DEBUG | 10 | Implemented |
| LOG_NONE | 0 | Implemented |
| TP | 0 | Implemented |
| HP | 1 | Implemented |
| SP | 2 | Implemented |
| TV | 3 | Implemented |
| UV | 4 | Implemented |
| SV | 5 | Implemented |
| NUM_REACTANTS | 0 | Implemented |
| NUM_PRODUCTS | 1 | Implemented |
| NUM_GAS | 2 | Implemented |
| NUM_CONDENSED | 3 | Implemented |
| NUM_ELEMENTS | 4 | Implemented |
| MAX_EQUATIONS | 5 | Implemented |
| TEMPERATURE | 0 | Implemented |
| PRESSURE | 1 | Implemented |
| VOLUME | 2 | Implemented |
| DENSITY | 3 | Implemented |
| ENTHALPY | 6 | Implemented |
| ENERGY | 7 | Implemented |
| ENTROPY | 8 | Implemented |
| GIBBS_ENERGY | 9 | Implemented |
| GAMMA_S | 10 | Implemented |
| FROZEN_CP | 11 | Implemented |
| FROZEN_CV | 12 | Implemented |
| EQUILIBRIUM_CP | 13 | Implemented |
| EQUILIBRIUM_CV | 14 | Implemented |
| VISCOSITY | 15 | Implemented |
| FROZEN_CONDUCTIVITY | 16 | Implemented |
| EQUILIBRIUM_CONDUCTIVITY | 17 | Implemented |
| FROZEN_PRANDTL | 18 | Implemented |
| EQUILIBRIUM_PRANDTL | 19 | Implemented |
| ROCKET_TEMPERATURE | 0 | Implemented |
| ROCKET_PRESSURE | 1 | Implemented |
| ROCKET_VOLUME | 2 | Implemented |
| ROCKET_DENSITY | 3 | Implemented |
| ROCKET_M | 4 | Implemented |
| ROCKET_MW | 5 | Implemented |
| ROCKET_ENTHALPY | 6 | Implemented |
| ROCKET_ENERGY | 7 | Implemented |
| ROCKET_ENTROPY | 8 | Implemented |
| ROCKET_GIBBS_ENERGY | 9 | Implemented |
| ROCKET_GAMMA_S | 10 | Implemented |
| ROCKET_FROZEN_CP | 11 | Implemented |
| ROCKET_FROZEN_CV | 12 | Implemented |
| ROCKET_EQUILIBRIUM_CP | 13 | Implemented |
| ROCKET_EQUILIBRIUM_CV | 14 | Implemented |
| MACH | 15 | Implemented |
| SONIC_VELOCITY | 16 | Implemented |
| AE_AT | 17 | Implemented |
| C_STAR | 18 | Implemented |
| COEFFICIENT_OF_THRUST | 19 | Implemented |
| ISP | 20 | Implemented |
| ISP_VACUUM | 21 | Implemented |
| ROCKET_VISCOSITY | 22 | Implemented |
| ROCKET_FROZEN_CONDUCTIVITY | 23 | Implemented |
| ROCKET_EQUILIBRIUM_CONDUCTIVITY | 24 | Implemented |
| ROCKET_FROZEN_PRANDTL | 25 | Implemented |
| ROCKET_EQUILIBRIUM_PRANDTL | 26 | Implemented |
| SHOCK_TEMPERATURE | 0 | Implemented |
| SHOCK_PRESSURE | 1 | Implemented |
| SHOCK_VELOCITY | 2 | Implemented |
| SHOCK_MACH | 3 | Implemented |
| SHOCK_SONIC_VELOCITY | 4 | Implemented |
| SHOCK_RHO12 | 5 | Implemented |
| SHOCK_RHO52 | 6 | Implemented |
| SHOCK_P21 | 7 | Implemented |
| SHOCK_P52 | 8 | Implemented |
| SHOCK_T21 | 9 | Implemented |
| SHOCK_T52 | 10 | Implemented |
| SHOCK_M21 | 11 | Implemented |
| SHOCK_M52 | 12 | Implemented |
| SHOCK_V2 | 13 | Implemented |
| SHOCK_U5_P_V2 | 14 | Implemented |
| SHOCK_VOLUME | 15 | Implemented |
| SHOCK_DENSITY | 16 | Implemented |
| SHOCK_M | 17 | Implemented |
| SHOCK_MW | 18 | Implemented |
| SHOCK_ENTHALPY | 19 | Implemented |
| SHOCK_ENERGY | 20 | Implemented |
| SHOCK_ENTROPY | 21 | Implemented |
| SHOCK_GIBBS_ENERGY | 22 | Implemented |
| SHOCK_GAMMA_S | 23 | Implemented |
| SHOCK_FROZEN_CP | 24 | Implemented |
| SHOCK_FROZEN_CV | 25 | Implemented |
| SHOCK_EQUILIBRIUM_CP | 26 | Implemented |
| SHOCK_EQUILIBRIUM_CV | 27 | Implemented |
| SHOCK_VISCOSITY | 28 | Implemented |
| SHOCK_FROZEN_CONDUCTIVITY | 29 | Implemented |
| SHOCK_EQUILIBRIUM_CONDUCTIVITY | 30 | Implemented |
| SHOCK_FROZEN_PRANDTL | 31 | Implemented |
| SHOCK_EQUILIBRIUM_PRANDTL | 32 | Implemented |
| DETONATION_P1 | 0 | Implemented |
| DETONATION_T1 | 1 | Implemented |
| DETONATION_H1 | 2 | Implemented |
| DETONATION_M1 | 3 | Implemented |
| DETONATION_GAMMA1 | 4 | Implemented |
| DETONATION_V_SONIC1 | 5 | Implemented |
| DETONATION_PRESSURE | 6 | Implemented |
| DETONATION_TEMPERATURE | 7 | Implemented |
| DETONATION_DENSITY | 8 | Implemented |
| DETONATION_ENTHALPY | 9 | Implemented |
| DETONATION_ENERGY | 10 | Implemented |
| DETONATION_GIBBS_ENERGY | 11 | Implemented |
| DETONATION_ENTROPY | 12 | Implemented |
| DETONATION_MACH | 13 | Implemented |
| DETONATION_VELOCITY | 14 | Implemented |
| DETONATION_SONIC_VELOCITY | 15 | Implemented |
| DETONATION_GAMMA | 16 | Implemented |
| DETONATION_P_P1 | 17 | Implemented |
| DETONATION_T_T1 | 18 | Implemented |
| DETONATION_M_M1 | 19 | Implemented |
| DETONATION_RHO_RHO1 | 20 | Implemented |
| DETONATION_FROZEN_CP | 21 | Implemented |
| DETONATION_FROZEN_CV | 22 | Implemented |
| DETONATION_EQUILIBRIUM_CP | 23 | Implemented |
| DETONATION_EQUILIBRIUM_CV | 24 | Implemented |
| DETONATION_M | 25 | Implemented |
| DETONATION_MW | 26 | Implemented |
| DETONATION_VISCOSITY | 27 | Implemented |
| DETONATION_FROZEN_CONDUCTIVITY | 28 | Implemented |
| DETONATION_EQUILIBRIUM_CONDUCTIVITY | 29 | Implemented |
| DETONATION_FROZEN_PRANDTL | 30 | Implemented |
| DETONATION_EQUILIBRIUM_PRANDTL | 31 | Implemented |
| R / constants.R | 8314.51 J/(kmol K) | Implemented |

## Unit helpers

Each takes one argument, value; no default; returns the same scalar/NumPy shape.

| Member | Conversion | Status |
|---|---|---|
| celsius_to_kelvin(value) | Convert °C to K. | Implemented |
| kelvin_to_celsius(value) | Convert K to °C. | Implemented |
| fahrenheit_to_kelvin(value) | Convert °F to K. | Implemented |
| kelvin_to_fahrenheit(value) | Convert K to °F. | Implemented |
| rankine_to_kelvin(value) | Convert °R to K. | Implemented |
| kelvin_to_rankine(value) | Convert K to °R. | Implemented |
| atm_to_bar(value) | Convert atmospheres to bar. | Implemented |
| bar_to_atm(value) | Convert bar to atmospheres. | Implemented |
| psi_to_bar(value) | Convert psi to bar. | Implemented |
| bar_to_psi(value) | Convert bar to psi. | Implemented |
| mmhg_to_bar(value) | Convert mmHg to bar. | Implemented |
| bar_to_mmhg(value) | Convert bar to mmHg. | Implemented |
| g_per_cm3_to_kg_per_m3(value) | Convert g/cm³ to kg/m³. | Implemented |
| kg_per_m3_to_g_per_cm3(value) | Convert kg/m³ to g/cm³. | Implemented |
| cm3_per_g_to_m3_per_kg(value) | Convert cm³/g to m³/kg. | Implemented |
| m3_per_kg_to_cm3_per_g(value) | Convert m³/kg to cm³/g. | Implemented |
| m_per_s_to_ft_per_s(value) | Convert m/s to ft/s. | Implemented |
| ft_per_s_to_m_per_s(value) | Convert ft/s to m/s. | Implemented |
| joule_to_cal(value) | Convert joules to calories. | Implemented |
| cal_to_joule(value) | Convert calories to joules. | Implemented |
| kj_per_kg_to_cal_per_g(value) | Convert kJ/kg to cal/g. | Implemented |
| cal_per_g_to_kj_per_kg(value) | Convert cal/g to kJ/kg. | Implemented |

## Other namespaces

cea.lib re-exports the same library API. The optional MATLAB module exposes
eq_solve, rocket_solve, shock_solve and detonation_solve, returning flattened
SimpleNamespace objects; its Matlab*Solution stub types carry the same fields
as their corresponding solutions. MATLAB wrappers, sample plotting and file
outputs are outside this adapter. Access raises a named NotImplementedError.
The native __all__ also accidentally exports os and importlib_resources;
these helper modules are not re-exported.

## Version, exports and optional wrappers

| Member | Status |
|---|---|
| __version__ | Differs: actual cea-keplr-rs distribution version |
| __all__ | Implemented supported API; excludes accidental native os/importlib_resources exports |
| constants | Implemented selector facts and R |
| units | Implemented all 22 helpers |
| lib / lib.libcea | Implemented re-exports; no Fortran shared library |
| matlab | NotImplementedError for its four wrappers and MATLAB-specific solution types |

MATLAB wrapper signatures are recorded for completeness. They are a separate
namespace, not root aliases. Their unsupported names raise NotImplementedError.

| Member | Signature | Status |
|---|---|---|
| eq_solve | `(eq_type: 'int', reactants: 'list[str]', T: 'Optional[float]' = None, H: 'Optional[float]' = None, S: 'Optional[float]' = None, U: 'Optional[float]' = None, P: 'Optional[float]' = None, V: 'Optional[float]' = None, T_reac: 'Optional[ScalarOrArray]' = None, fuel_amounts: 'Optional[VectorLike]' = None, oxid_amounts: 'Optional[VectorLike]' = None, moles: 'bool' = False, of_ratio: 'Optional[ScalarOrArray]' = None, phi: 'Optional[ScalarOrArray]' = None, r_eq: 'Optional[ScalarOrArray]' = None, pct_fuel: 'Optional[ScalarOrArray]' = None, only: 'Optional[list[str]]' = None, omit: 'Optional[list[str]]' = None, insert: 'Optional[list[str]]' = None, trace: 'Optional[float]' = None, transport: 'bool' = False, ions: 'bool' = False) -> 'SimpleNamespace'` | NotImplementedError: MATLAB wrapper |
| rocket_solve | `(reactants: 'list[str]', pc: 'SupportsFloat', pi_p: 'Optional[ScalarOrArray]' = None, subar: 'Optional[ScalarOrArray]' = None, supar: 'Optional[ScalarOrArray]' = None, T_reac: 'Optional[ScalarOrArray]' = None, fuel_amounts: 'Optional[VectorLike]' = None, oxid_amounts: 'Optional[VectorLike]' = None, moles: 'bool' = False, of_ratio: 'Optional[SupportsFloat \| VectorLike]' = None, phi: 'Optional[SupportsFloat \| VectorLike]' = None, r_eq: 'Optional[SupportsFloat \| VectorLike]' = None, pct_fuel: 'Optional[SupportsFloat \| VectorLike]' = None, only: 'Optional[list[str]]' = None, omit: 'Optional[list[str]]' = None, insert: 'Optional[list[str]]' = None, trace: 'Optional[float]' = None, transport: 'bool' = False, ions: 'bool' = False, iac: 'bool' = True, n_frz: 'Optional[int]' = None, hc: 'Optional[SupportsFloat]' = None, tc: 'Optional[SupportsFloat]' = None, mdot: 'Optional[SupportsFloat]' = None, ac_at: 'Optional[SupportsFloat]' = None, tc_est: 'Optional[SupportsFloat]' = None) -> 'SimpleNamespace'` | NotImplementedError: MATLAB wrapper |
| shock_solve | `(reactants: 'list[str]', T0: 'SupportsFloat', p0: 'SupportsFloat', *, u1: 'Optional[SupportsFloat \| VectorLike]' = None, Mach1: 'Optional[SupportsFloat \| VectorLike]' = None, fuel_amounts: 'Optional[VectorLike]' = None, oxid_amounts: 'Optional[VectorLike]' = None, moles: 'bool' = False, of_ratio: 'Optional[SupportsFloat \| VectorLike]' = None, phi: 'Optional[SupportsFloat \| VectorLike]' = None, r_eq: 'Optional[SupportsFloat \| VectorLike]' = None, pct_fuel: 'Optional[SupportsFloat \| VectorLike]' = None, only: 'Optional[list[str]]' = None, omit: 'Optional[list[str]]' = None, insert: 'Optional[list[str]]' = None, trace: 'Optional[float]' = None, transport: 'bool' = False, ions: 'bool' = False, reflected: 'bool' = True, incident_frozen: 'bool' = False, reflected_frozen: 'bool' = False) -> 'SimpleNamespace'` | NotImplementedError: MATLAB wrapper |
| detonation_solve | `(reactants: 'list[str]', T1: 'SupportsFloat', p1: 'SupportsFloat', *, fuel_amounts: 'Optional[VectorLike]' = None, oxid_amounts: 'Optional[VectorLike]' = None, moles: 'bool' = False, of_ratio: 'Optional[SupportsFloat \| VectorLike]' = None, phi: 'Optional[SupportsFloat \| VectorLike]' = None, r_eq: 'Optional[SupportsFloat \| VectorLike]' = None, pct_fuel: 'Optional[SupportsFloat \| VectorLike]' = None, only: 'Optional[list[str]]' = None, omit: 'Optional[list[str]]' = None, insert: 'Optional[list[str]]' = None, trace: 'Optional[float]' = None, transport: 'bool' = False, ions: 'bool' = False, frozen: 'bool' = False) -> 'SimpleNamespace'` | NotImplementedError: MATLAB wrapper |

## Limits and bounded exceptions

The adapter changes units and data layout, not numerical methods. The original
cea_rs API and solver validation remain in force. There is no fallback to NASA.
Before the first successful solve, or after a failed solve, reading solution
properties raises RuntimeError; converged remains False. In contrast, native
CEA initializes numeric workspaces and can retain an older failed solution.
Its printed report, CLI, sample plotting and file outputs are not reproduced.

- RP-1311 example 12, equilibrium area 75 only: CO/OH and H2 mole fractions
  exceed 1e-4 relative because native stops at area 74.99733854. The tests bound
  these exact species/area cells at 1.8e-4, including corresponding mass fractions
  and nj. H2 is newly visible in this API comparison because the historical
  example-12 report uses mass fractions below the reporting threshold for H2.
- RP-1311 example 7 at 1000 m/s: native returns LAST_VALID_SOLUTION=9 and
  converged=False. Rust raises ConvergenceError, leaving no readable stale state.
- Reflected-frozen station 5: native Cp/Cv fields are uninitialized. Rust uses
  current frozen composition/temperature for Cp and Cv, checked independently.
  With transport enabled, native cp_fr is populated and compares normally.
  Native gamma and sound speed retain the station-2 Cp basis in both backends.
  Shock cv_eq is zero at initial/frozen stations, following the native sentinel.
  Transport is zero at station 1; equilibrium conductivity/Prandtl are zero
  at frozen shock stations, matching NASA.
- Frozen rockets: native Cv can retain an earlier value or become zero. Rust
  returns current Cv; tests check Cp-R*n independently and bound native values
  to zero or an earlier equilibrium Cv. The chamber equilibrium Cp is preserved
  for n_frz=1, including NASA subsonic schedule behavior.
- RP-1311 example 9 at n_frz=3: native error8 retains seven stations. Rust
  returns ten complete verified stations. Tests require that exact native failure
  and independently check Rust stagnation enthalpy. Other tested FAC freeze
  indices (4/5) compare at 1e-4.
- Air/Jet-A(L), O/F25, 20 MPa, area80: NASA reports success with 14–15 kJ/kg
  energy error and 18–19% mass-flux error using its own reported properties.
  Rust reaches water coexistence at 273.121944 K and 9971.217844 Pa, with
  independently checked entropy/energy/mass closure. This exact case is a
  bounded exception, not a comparison against NASA invalid output.
- Condensed transitions: native error8 at Ar/Hg HP, Ar/water TP273.15 K, and
  Ar/NaCN TP288.5 K is bounded by exact input and Rust temperature/phase values
  from independent checks in condensed-phase-insertion.md. Ar/Br2 TP262 K
  omits the solid in native because of its reversed interval; Rust returns
  solid mass fraction 0.943263839016. RP14 water condensation and RP13 BeO
  coexistence compare normally. Cv_eq=0 at Cp_eq=0 phase coexistence follows
  the native reporting convention, avoiding a singular derivative query.
- NASA ShockSolution.n/nj/ln_nj call a missing _get_size method and raise
  AttributeError. Rust returns n/nj with the documented shapes, checked against
  native fraction/molar-mass sums. Untruncated frozen ln_nj is unavailable here
  and explicitly raises NotImplementedError; log(nj) is not substituted.
  Element-reduced pools also raise if an untruncated log is unavailable.
- Transport Cp is the native gas-only value; unsupported fits/transport failures
  still raise this port’s typed errors. Native estimates can sometimes continue.
- FD derivatives of an assigned H/U/S with respect to the other input or
  weights are mathematically zero. At h=1e-5 both native and Rust amplify
  forward-solve residuals: tests bound these exact HP/SP/UV/SV responses at
  1e-3 (kJ-based H/U) or 1e-5 (kJ-based S), and compare other FD responses
  normally at 1e-4 relative. Analytic responses agree directly.
- Other unsupported features raise NotImplementedError by name: arbitrary
  initial guesses, smooth truncation, fitted named-reactant overrides (which
  would require preserving/altering their polynomial fits), n_frz=-1 frozen
  reactant rockets, frozen detonation, custom binary databases, MATLAB wrappers
  and verbose derivative diagnostics. Analytic derivative limits of the core
  (ions, element-reduced systems and UV/SV with condensates) still apply.
- Nominal fit metadata do not extend this port’s validated domains. In
  particular N2O at 298.15 K can still be rejected where native extrapolates.
  The complete pre-existing limits remain in implementation-matrix.md.

Comparison policy: relative tolerance 1e-4 on thermodynamic/performance outputs
and species with mole fraction >=1e-5. Exact zeros use explicit, small absolute
allowances. Derivative near-zero noise has separate absolute bounds. Exceptions
are restricted to named inputs/fields with independent identities or exact
failure codes; no test is skipped. All helper calls are identical for both
backends; the adapter never imports the native package.

## Validation on 2026-11-17

macOS arm64; Rust 1.99.1, maturin 1.16.0, NumPy 2.5.5. The installed wheel
was tested outside the source tree with isolated Python imports.

| Problem type | Worst relative difference | Member |
|---|---:|---|
| TP | 3.7932151e-10 | cv_eq |
| HP | 3.560231e-11 | cv_eq |
| SP | 3.6130657e-11 | cv_eq |
| TV | 2.1268467e-09 | cv_fr |
| UV | 3.5188364e-11 | cv_eq |
| SV | 3.5204482e-11 | cv_eq |
| IAC equilibrium | 7.0467525e-05 | mass_fractions.NO |
| IAC frozen | 1.5181817e-05 | P |
| FAC equilibrium | 6.8592549e-05 | mass_fractions.H |
| FAC frozen | 3.1813676e-05 | P |
| shock | 1.8963376e-06 | cv_fr |
| detonation | 4.1122016e-05 | energy |
| condensed TP | 2.0456017e-09 | nj |
| mixture | 7.5634463e-16 | calc_property 7 ndarray |
| derivatives | 2.2131531e-12 | TV.dnj_dstate2 |
| native derivative FD | 3.7132636e-07 | UV.dH_dstate2_fd |

The separate RP12 exception is 1.73451869e-4 relative (CO mass fraction at
area75), below its named 1.8e-4 bound. All other table rows meet 1e-4.
Fixed-constraint FD noise and invalid native states are bounded separately
as described above; the table does not claim parity for those values.

| Check | Result |
|---|---|
| cargo test --release | 70 passed, 0 failed, 1 ignored (existing diagnostic) |
| Clippy, default, --release --all-targets -- -D warnings | Passed |
| Clippy, --features python, same options | Passed |
| cargo fmt --check | Passed |
| uvx maturin build --release | Passed; cp312-abi3 macosx_11_0_arm64 |
| Wheel contents | Python adapter, stubs, py.typed and optional nasa/NumPy extra verified |
| Python 3.12.14, neither NumPy nor cea installed | test_wheel: 8 passed; 0 failures/skips |
| Python 3.15.1, neither NumPy nor cea installed | test_wheel: 8 passed; 0 failures/skips |
| Python 3.12.14, with cea==3.3.4 and NumPy | test_nasa_api: 9 passed; 0 failures/skips |
| Python 3.15.1, with NumPy, cea confirmed absent | test_nasa_api: 9 passed; 0 failures/skips |

[Raw results, named exceptions and exact check commands](nasa-python-api-results.json).
The fallback suite checks 41 stored NASA-only snapshots, covering all successful
compared rocket schedules, shock flag combinations, detonation cases and RP14
condensation points.
Known stale/uninitialized native heat capacities are excluded from snapshots;
their independent identities and live-native exception checks remain in place.
[Snapshot provenance](../tests/data/nasa/README.md) records the source and inputs.
CI installs the base wheel with `--no-deps`, runs its tests, then installs the
`nasa` extra into a separate directory for adapter tests without NASA itself
on Python 3.12 and 3.15.
Linux/Windows wheels and production consumers were not tested locally.
