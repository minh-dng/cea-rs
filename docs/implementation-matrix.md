# Implemented scope and remaining work

Version 0.1.0a1 is an experimental, separately named distribution. Full NASA CEA
parity, Python API compatibility and consumer replacement are not established.

| Capability | Rust status | Evidence or remaining gate |
|---|---|---|
| NASA nine-coefficient Cp, H, S | Implemented | Named product/reactant fits and reference-only reactants; exact NASA data identity |
| Custom formulas and enthalpy | Implemented | Custom N2O/ABS three-component fuel comparisons; caller supplies consistent formula/molar mass |
| TP, HP, SP, TV, UV, SV | Implemented | Both H2/O2 and CH4/O2/N2 compared with native 3.3.4 |
| Condensed products | Experimental active set | Graphite comparison; general phase coexistence and metal systems unproven; 11 invalid data records omitted |
| Ions | Experimental neutral balance | H/H2/H+/e- declared pool at 10000 K; broad ionized mixtures unproven |
| IAC equilibrium nozzle | Implemented | Chamber, sonic throat, one supersonic exit; five propellant sets at two operating points |
| Frozen at chamber / throat | Implemented | Native n_frz=1/2 comparisons; fixed species and phases; equilibrium Cp unavailable for frozen composition |
| Finite-area chamber (FAC) | Unsupported | Conservation closure, chamber area/mass-flux API and native regression cases required |
| Incident/reflected shocks | Unsupported | Equilibrium/frozen Rankine-Hugoniot solver and native examples required |
| Chapman-Jouguet detonation | Unsupported | CJ sonic closure, initial-state model and native examples required |
| Transport | Unsupported | NASA transport fits, missing-species rules, gas/condensed mixing and equilibrium transport required |
| Arbitrary exit pressures/subsonic stations | Unsupported | Currently one area-specified supersonic exit only |
| Input-deck CLI, C/Fortran/Excel/MATLAB bindings | Unsupported | No legacy input or ABI compatibility layer |
| NASA Python classes and options | Unsupported | No Mixture/EqSolution/RocketSolution mutation API, equivalence-ratio utilities, trace/smoothing controls or EqDerivatives API |
| Analytic sensitivities | Unsupported | Internal equilibrium heat capacity/gamma do not imply a public derivative API |
| Ambient/separated nozzle thrust | Unsupported | Outputs are ideal vacuum velocity; no consumer thrust-efficiency model is inherited |

## Units and failures

All solver inputs and state outputs use SI: K, Pa, kg/m³, J/kg, J/(kg K), and
m³/kg. Molar mass is numerically g/mol = kg/kmol. Reactant `thermo` uses J/mol
and J/(mol K); custom enthalpy is J/mol at its declared reference temperature.
`R` is NASA CEA's 8314.51 J/(kmol K), deliberately retained for comparison.
Rocket velocities and c-star are m/s. Divide vacuum velocity by 9.80665 for
specific impulse in seconds. Chamber area ratio is None for an IAC model.

Invalid weights, formulas, product pools and nonfinite inputs raise ValueError.
DomainError rejects out-of-fit reactant enthalpy, appreciable out-of-fit
products, and solved temperatures outside 300–20000 K. A reference-only
reactant has no heat capacity and is accepted only at its exact reference
state, with 1e-7 K rounding tolerance. No low-temperature fit is invented.

NASA N2O gas starts at 300 K. Hybrid's canonical 298.15 K request therefore
raises DomainError here; native CEA extrapolates. The N2O/ABS comparisons use
300 K N2O and 298.15 K custom fuel references. This difference alone prevents
silent use as a Hybrid substitute. Cold-flow Hybrid simulations are not
validated by these thermochemical tests.

The attributed Newton implementation may evaluate fits outside their ranges
while iterating. Returned species above 1e-11 mole fraction must be within
range. Any nonzero out-of-fit trace species below that cutoff are listed in
`State.trace_extrapolated_species`; their extrapolated trace contribution is
a known limit, not an extended scientific validity claim.

ConvergenceError rejects exhausted iteration budgets, singular thermodynamic
derivatives, nonfinite/nonphysical output, excessive element/charge/constraint
residuals and failed nozzle sonic/mass/energy closure. No stale result or
native fallback is returned. Energy-constraint residuals use
`max(abs(target), RT/M, 1 J/kg)` so near-zero reference enthalpy remains
meaningful. Other constraints use `max(abs(target), 1 SI unit)`.

Mole fractions follow NASA's reporting convention: their denominator includes
condensed moles. They are not gas-only mixing fractions. `M=1/n_gas` is the
mass-based effective molar mass. Python equilibrium Cp is None in frozen
states. The Rust field is valid only when `State.equilibrium` is true.

## Production and performance gate

Keep native CEA as the production default. The release benchmark shows the
current Rust implementation is slower, especially the bracketed nozzle's
many independent equilibrium solves. Optimization requires measured profiles,
reused immutable species/model data and better root algorithms, with the same
native and residual tests retained. Language alone does not establish speed.
