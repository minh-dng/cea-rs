# Implemented scope and remaining work

Version 0.1.0a3 is an experimental, separately named distribution. Full NASA CEA
parity and consumer replacement are not established. The supported NASA Python
API subset is described in [the member inventory](nasa-python-api.md).

| Capability | Rust status | Evidence or remaining gate |
|---|---|---|
| NASA nine-coefficient Cp, H, S | Implemented | NASA gas extensions, reactant reference/extrapolation rules and exact source-data identity; [parity evidence](equilibrium-parity.md) |
| Custom formulas and enthalpy | Implemented | Custom N2O/ABS three-component fuel comparisons; caller supplies consistent formula/molar mass |
| TP, HP, SP, TV, UV, SV | Implemented | NASA 160–22000 K bounds; RP-1311 examples 1–5 and 14 plus native regressions; extrapolated products exposed on State |
| Condensed products | Experimental active set | Compiled phase numbering/coexistence and all 11 previously skipped records restored; aluminized and low-temperature native comparisons pass; broader metal systems unproven. Cold HP/SP gas-only failures retry from a converged 160 K condensed state. Inverted condensed intervals retain their upper phase limit without using the invalid lower bound. Fixed-T transitions select the lower-Gibbs sibling without a singular coexistence pair; Gibbs energies are compared only inside the candidate phase's own fit range, because an extrapolated fit can look falsely stable (NaCN(III) below 288.5 K); variable-T coexistence is preserved. Failed fixed-T solves retry once with nonzero trace pivots, in-range inserted phases and a gas-untruncation retry before singular phase removal. This resolves NaCN at 300 K, insert NaCN(II) at 300/835 K, and private ice-to-liquid TP continuation while preserving every previously successful result. The Hg cold-start, reversed Br2 interval, exact-transition and fixed-T recovery cases have first-principles regressions; see [condensed-phase checks](condensed-phase-insertion.md). |
| Ions | Experimental neutral balance | H/H2/H+/e- declared pool at 10000 K; broad ionized mixtures unproven |
| IAC equilibrium nozzle | Implemented station schedules | Private temperature/composition/phase continuation; RP-1311 13 with its trace=1e-10 now passes all six stations at 1e-5 relative, alongside 8 and 12. Air/Jet-A(L), O/F=25, 0.7 MPa and area=200 passes at 1e-5 using a refined native sonic throat. At 20 MPa and area=80, corrected phase adjacency reaches a verified ice/liquid exit at 273.121944 K and 9971.217844 Pa. An independent native 3.3.4 rerun reports success but fails entropy, energy, mass-flux and water-saturation checks; its standalone SP solve fails with error 8. All 89 neighboring cases retain identical results. See [measured evidence](nozzle-continuation.md#20-mpa-water-phase-audit) and [rocket evidence](rocket-parity/README.md) |
| Frozen at chamber / throat | Implemented station-index freezing | Native n_frz=1/2/3/5 comparisons, schedule omissions and gas/condensed fit stops; fixed species and phases; equilibrium Cp unavailable for frozen composition |
| Finite-area chamber (FAC) | Experimental | ac_at or mass flux per chamber area; injector/infinity/combustor/throat stations; native RP-1311 9/10 pass; early-freeze and sonic-chamber differences recorded in [rocket evidence](rocket-parity/README.md) |
| Incident/reflected shocks | Implemented with fit-domain limits | SI Rust/Python API; all four equilibrium/frozen combinations; independent mass/momentum/energy closure at 5e-5. Native 3.3.4 goldens: RP-1311 example 7 at both pressures/all seven speeds (1000 m/s equilibrium fails in native and returns an error here), four reflected H2/O2/Ar cases, N2/O2 air at Mach 3/6, frozen/frozen air at Mach 10, and all four Mach-10 modes with an explicit N2/O2/NO/N/O pool; worst relative difference 7.93e-6 against 1e-4 tolerance. Unrestricted air Mach-10 equilibrium reflection needs out-of-fit N2O; trace-CO2 air also exposes existing equilibrium residual failures. Frozen Cp is evaluated at the station temperature; reflected-frozen sound speed/gamma preserve NASA's incident-Cp convention. Transport excluded. |
| Chapman-Jouguet detonation | Implemented equilibrium CJ, SI | Fixed-composition unburned gas fits; NASA's damped iteration; fresh burned state with mass/momentum/energy/sonic checks. Native 3.3.4: RP-1311 example 6, CH4/O2/N2 air and acetylene/O2 at phi=0.5/1/2 meet 1e-4 relative. An extra named-Air case needs 3e-4 for near-zero internal energy; shared equilibrium parity remains open. Frozen detonation is unsupported in native CEA and here. Native Python rounds T1 to float32; goldens use its effective T1, with a separate exact 298.15 K fit check. |
| Transport | Implemented for State and rocket stations in Rust/Python | NASA gas selection, missing-data estimates, ion interactions, frozen/equilibrium Cp, conductivity and Prandtl; condensed species excluded, so transport Cp is gas-only. Native 3.3.4 comparisons include restricted gas pools and pass 1e-4 relative. Rocket `transport=True` (Rust `RocketOptions::transport`) fills `Station.transport` at every station from that station's composition; frozen stations use the frozen composition, as in NASA. All seven RP-1311 example 11 stations match NASA's reference output for all seven transport rows (worst 2.47e-5), and frozen n_frz=1/2 throat and exit stations match native. Unlike NASA, the station's thermodynamic Cp is not replaced by the transport Cp. Also unlike NASA, which fills estimates and continues, a transport failure at any station fails the whole rocket call. Finite-area-chamber stations and frozen stations with condensed species have no stored native golden. |
| Arbitrary exit pressures/subsonic stations | Implemented | Any number of pi_p/subar/supar stations; native pressure and area Newton iterations; large-subar native stopping defect recorded in [rocket evidence](rocket-parity/README.md) |
| NASA input decks | Implemented for sample keywords | Rust typed parser, case runner and Python `run_deck`; arbitrary dataset order, mixed units, custom formulas/enthalpies, mixture and station schedules; unknown keywords fail by name. All 14 decks run; strict comparison results below |
| Input-deck CLI, C/Fortran/Excel/MATLAB bindings | Unsupported | Library/Python deck entry points are available; no CLI, report formatter or ABI compatibility layer |
| NASA Python classes and options | Supported subset in `cea_rs.nasa` | Mixture, Reactant, mutable equilibrium/rocket/shock/detonation solutions and solvers, transport, sensitivities, root eq_solve, constants and units. Preserves NASA units, product order, float32 scalar boundaries and NumPy shapes. Smoothing, arbitrary initial guesses, binary databases, frozen-reactant rockets, untruncated frozen logs and MATLAB wrappers raise named NotImplementedError. [Complete member status and bounded differences](nasa-python-api.md). No solver numerical changes |
| Analytic sensitivities | Converged hard truncation | Solver.derivatives returns SI derivatives of T, gas n, species amounts, H/U/G/S by both state inputs and raw weights. TP/HP/SP native goldens for H2/O2 and CH4/O2/N2 at 2500 K, 1 MPa; worst relative difference 6.98e-7 (tolerance 1e-4). All six gas modes cross-checked; finite differences reject changed product pools/phases. Core cutoff and original gas logs are retained privately on State, including the residual cutoff before the trace-reporting pass. AP/Al at 300 K, 1 bar takes singular recovery to cutoff 80; all bulk sensitivities agree with central differences at 1e-5 relative and species at 1e-4 plus 1e-12 kmol/kg per input unit for trace noise. Analytic ions, element-reduced systems and UV/SV with condensates return errors; smoothing is unsupported |
| Ambient/separated nozzle thrust | Ideal ambient-matched performance | Per-station Cf, Isp and Isp_vacuum in SI; no assigned ambient pressure, separation or consumer thrust-efficiency model |
| Extra state properties | Implemented | Gibbs energy, frozen/equilibrium Cv, specific volume and NASA M/MW; native equilibrium and graphite goldens at 1e-4 relative |

## NASA sample decks

Every main-interface deck and both sample inputs at NASA commit `2f79a647`
are in `tests/data/nasa`, with unchanged original reports and NASA's notice.
Tests parse report values directly. Native 3.3.4 reports from an independently
rebuilt pinned executable provide an additional oracle. All 55 main-interface
cases solve. Nonzero values use `1e-4` relative tolerance, with composition
>=1e-5; zero unburned reference enthalpy uses 1e-3 J/kg absolute tolerance.
Example 6 uses the complete historical report because the new native printer
mislabels its unburned table. No solver numerics were changed.

| Deck | Primary comparison | Worst relative difference | Historical reference issue |
|---|---|---:|---|
| 1 | Pass | 1.61392e-9 | Two first-case carbon fractions differ by 1.10248e-4 |
| 2 | Pass | 6.74146e-8 | Pass |
| 3 | Pass | 4.40289e-10 | Pass |
| 4 | Pass | 1.39425e-10 | Pass |
| 5 | Pass | 1.48335e-5 | Pass |
| 6 | Pass, historical long report | 4.11214e-5 | Modern unburned/burned print defect |
| 7 | Pass | 8.81075e-6 | Old mmHg conversion; worst 1.68507e-4 |
| 8 | Pass | 6.96275e-5 | Subsonic performance; worst 1.14411e-4 |
| 9 | Pass | 6.85920e-5 | Combustor performance; worst 1.14472e-4 |
| 10 | Pass | 6.68535e-5 | Combustor area/performance; worst 1.39245e-4 |
| 11 | Pass | 3.20473e-5 | Subsonic performance; worst 1.00125e-4 |
| 12 | **Fail**, equilibrium CO/OH at area 75 | 1.73452e-4 | Original report contains frozen-only input; frozen comparison passes |
| 13 | Pass | 2.12241e-5 | Pass |
| 14 | Pass | 2.04574e-9 | Pass |
| Standalone sample 1 | Pass | 1.61392e-9 | Equivalent to example 1 |
| Combined RP-1311 | 56/57 solve; known 1000 m/s shock failure; two example 12 trace differences | 1.73455e-4 among successes | All other successful cases compared directly |
| Frozen detonation | Unsupported | — | Explicit error, as in native 3.3.4 |

Example 12 is a nozzle stopping difference: NASA reports area 74.99733854
for requested 75, while Rust reaches 75.00000000. At NASA's exact T/P the
public TP solver matches CO/OH within 6.4e-9 relative. The combined sample's
1000 m/s equilibrium shock fails in NASA as well; its printed retained values
are not accepted. The regression suite bounds explicitly named known cells
and rejects any new discrepancy. Green tests do not imply full deck parity.
[Deck API, comparison policy, exact defects and existing-test overlap](input-decks.md).

## Units and failures

The Rust and original `cea_rs` APIs use SI: K, Pa, kg/m³, J/kg, J/(kg K), and
m³/kg. Molar mass is numerically g/mol = kg/kmol. Reactant `thermo` uses J/mol
and J/(mol K); custom enthalpy is J/mol at its declared reference temperature.
The separate `cea_rs.nasa` adapter uses NASA's bar/kJ conventions instead.
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

Consumers configure their production backend. The 0.1.0a2 benchmark found
Rust slower, especially the bracketed nozzle's many independent equilibrium
solves; 0.1.0a3 has not been benchmarked. Optimization requires measured profiles,
reused immutable species/model data and better root algorithms, with the same
native and residual tests retained. Language alone does not establish speed.

## Tracked continuation work

- [Measure and improve Rust solve speed before consumer adoption](https://gitlab.com/keplr-space/cea-rs/-/work_items/1)
- [Implement finite-area chamber and additional nozzle station controls](https://gitlab.com/keplr-space/cea-rs/-/work_items/2)
- [Implement incident and reflected shock solvers](https://gitlab.com/keplr-space/cea-rs/-/work_items/3)
- [Implement Chapman-Jouguet detonation solver](https://gitlab.com/keplr-space/cea-rs/-/work_items/4)
- [Implement NASA transport properties and mixing rules](https://gitlab.com/keplr-space/cea-rs/-/work_items/5)
- [Plan explicit NASA API and input-deck compatibility](https://gitlab.com/keplr-space/cea-rs/-/work_items/6)
- [Broaden condensed-phase and ion numerical acceptance](https://gitlab.com/keplr-space/cea-rs/-/work_items/7)
