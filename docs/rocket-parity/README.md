# Rocket port validation

The rocket port extends the existing SI API with pressure, subsonic-area and
supersonic-area station lists; finite-area chambers specified by contraction
ratio or mass flow per chamber area; station-index freezing; and per-station
sonic velocity, c-star, thrust coefficient, ambient-matched Isp and vacuum Isp.
Assigned stagnation enthalpy `hc` is J/kg, assigned temperature `tc` and its
optional estimate `tc_est` are K. Isp values are velocities in m/s, as in NASA's
SI output, rather than seconds. Station transport is covered in
[the implementation matrix](../implementation-matrix.md).

```python
from cea_rs import Solver
solver = Solver(["H2(L)", "O2(L)"])
r = solver.rocket([1, 5.55157], [20.27, 90.17], 53.3172e5,
                  iac=False, ac_at=1.58,
                  pi_p=[10, 100, 1000], supar=[25, 50, 75])
for station in r.stations:
    print(station.name, station.state.temperature_k, station.isp_vacuum_m_s)
```

`Solver::rocket` and the existing positional Python signature still return the
chamber, throat and one exit. Rust's extended entry point is
`Solver::rocket_stations(..., &RocketOptions { ... })`. The result type is
unchanged and gains `stations`, `throat_index`, `complete` and `stop_reason`.
`throat_index` is zero-based; `n_frz` is NASA's one-based station index.

IAC station order is infinity, throat, then exits grouped by pi_p, subar, supar.
FAC adds injector, infinity and combustor end before the throat. `chamber` means
the FAC combustor end, and the IAC infinity station. The auxiliary FAC infinity
station has sonic output zero, as NASA reports; its State still has a physical
sound speed. With no exit requests, `exit` denotes the throat.

Freezing uses fixed amounts and phases from the requested station. NASA's
omission rules are retained: later pressure points above the freeze pressure,
later supersonic areas no larger than the freeze area when n_frz>=3, and subar
schedules when n_frz>1 are omitted. Invalid indices are errors, rather than
silently disabling frozen flow. Negative n_frz (NASA's separate reactant-frozen
mode) is not a station index and is not implemented.

Frozen gas calculations stop below 0.8 times the minimum lower gas fit bound
among species with amount above 1e-12 kmol/kg. Condensed phases stay within the
reference polynomial interval plus NASA's 50 K guard band. This removes the
old 300 K rocket floor. NASA's small extrapolation allowance is explicit:
`trace_extrapolated_species` lists every out-of-fit frozen species, including
appreciable ones in the allowed guard band. A station-list solve returns only
completed, verified stations and `complete=False` on a fit-limit stop. Legacy
single-exit calls raise DomainError if they cannot reach their requested exit.

## Native comparisons

Oracle: separately installed native `cea==3.3.4` from PyPI, with numpy, on macOS
arm64 and Python 3.15.0. NASA source inspected read-only at
`2f79a647e85737742fb08a127da00d1ea1d069da`. No NASA source-repository scripts were
executed and no native result enters the Rust solver. Run date: 2026-11-15.
Full station values, input models and differences are in [comparison.json](comparison.json).

All 15 requested state/flow/performance fields use strict relative error 1e-4,
with absolute error at zero. Rust's committed regression goldens add an
absolute floor of 1e-9 solely at zero-valued outputs. The native goldens were
produced by `scripts/compare_rocket_parity.py` and stored verbatim as decimal
values in `tests/rocket-native.txt`. The script exits nonzero while any listed
case fails; its known failures have not been suppressed.

| RP-1311 case | Result | Worst relative difference |
|---|---|---:|
| 8, IAC pressure/subsonic/supersonic stations | Pass | 2.25104e-5 |
| 9, FAC ac_at=1.58 | Pass | 2.20796e-5 |
| 10, FAC mdot=1333.9 kg/(s m²) | Pass | 2.44668e-5 |
| 11, named Li(cr) input | Pass | 9.42705e-6 |
| 11, explicit identical Li reference | Pass | 9.42702e-6 |
| 12, equilibrium | Pass | 4.60226e-5 |
| 12, frozen at throat | Pass | 1.52142e-5 |
| 13, BeO condensed phases, trace=1e-10 | Pass | 1.89622e-6 |

Twenty of the 25 comparison cases pass. These include freeze points 1,
2, 3 and 5, FAC freeze points 4 and 5, omitted/reordered frozen schedule points,
assigned hc/tc/tc_est, near-stagnation pressure exits, subar during chamber
freezing, and gas/condensed fit-limit partial stops. The existing independent
146-row comparison still passes, with worst scaled difference 5.83984e-5;
its metadata and summary are retained in
[existing-comparison-summary.json](existing-comparison-summary.json).

Every failed or unverified case is listed here:

- **fac-fr-1:** example 9 inputs, supar=[5,25,50], n_frz=1. NASA computes the
  combustor with the equilibrium throat, then replaces the throat with a
  frozen throat without updating chamber closure. Its combustor/throat mass
  flux ratio is 0.9067058557. Rust closes the frozen FAC throughout, keeping
  ac_at and injector momentum consistent. The worst native difference is
  8.47504e-2; this case is not accepted as native parity.
- **fac-fr-2 and fac-fr-3:** the same FAC case with n_frz=2 or 3. Native CEA
  reports CEA_NOT_CONVERGED and returns nonfinite throat/performance fields.
  Rust returns finite states satisfying sonic, area, mass, energy and injector
  momentum checks. There is no valid native result against which to verify
  these two cases.
- **subar-range:** example 8 reactants and pressure, subar=[1.001,1.01,2,10,100].
  Native's small-pressure-step stopping rule returns Ae/At=210.6552002798 for
  the last requested ratio 100. Rust bounds the geometric error before
  accepting that stopping rule and returns 100.0011376867. The largest relative
  difference is 1.10656. This case does not match native; conservation and
  assigned-geometry regression checks cover the Rust result.
- **fac-sonic:** example 9 with ac_at=1 and supar=[5,25]. Native's combustor
  Mach is 0.9965482473; Rust closes the coincident chamber/throat at Mach 1.
  The largest native difference is 3.68836e-3 in combustor pressure. This
  limiting geometry is checked by conservation, not accepted as native parity.

The guards keep tighter sonic and supersonic-area closure than native's
stopping tolerances. Subsonic calculations retain the native pressure-step
criterion only when the geometric error is below 5e-4. All returned moving
stations must close mass and energy within 2e-6; throat Mach must close within
2e-6. FAC injector momentum is solved to 2e-8. The edge-case tests enforce
geometry and conservation independently of native outputs.

## Timing

Release wheel built with Rust 1.99.1. Each call is the existing three-station
IAC solve at pc=3 MPa, Ae/At=10, equilibrium flow. Baseline is commit
64cc1c6 before the rocket edits. Construction is excluded. Native allocates a
fresh solution per call. The repeated figure is the median of 50 separately
timed calls after five warmups; the first-call columns retain a single measured
solve before warmup. Times are milliseconds.

| Propellant | Native first before / after | Rust first before / after | Native median after | Rust median before / after | Median speedup |
|---|---:|---:|---:|---:|---:|
| H2/O2 | 0.098333 / 0.082166 | 1.206708 / 0.418667 | 0.054583 | 1.157229 / 0.385792 | 3.00× |
| LOX/RP-1 | 0.936334 / 0.914750 | 9.556250 / 2.971084 | 0.836041 | 9.149687 / 2.888438 | 3.17× |
| N2O4/Aerozine-50 | 1.280958 / 1.307792 | 12.153250 / 4.297333 | 1.136583 | 11.603979 / 4.206500 | 2.76× |

Raw measurements are in [benchmark-before.json](benchmark-before.json) and
[benchmark-after.json](benchmark-after.json). The new pressure iteration makes
Rust about 2.76–3.17 times faster in these equilibrium cases, while it remains
3.45–7.07 times slower than native. CPU scheduling, thermal state and OS cache
were not controlled. These measurements do not establish a native speed win.

## Checks and integration hooks

Passed: cargo +1.99.1 test --locked (28 tests, one pre-existing ignored
diagnostic), cargo +1.99.1 clippy --features python -- -D warnings,
installed-wheel tests (four), the 146-row comparison, and Ruff on the changed
Python wrapper and new comparison script. The rocket comparison intentionally
fails overall because the seven cases above are still failed or unverified.

Reproduction requires an external throwaway environment with cea==3.3.4,
numpy and the installed release wheel:

```sh
python scripts/compare_rocket_parity.py --output /tmp/rocket-comparison.json
python scripts/compare_native.py --installed /path/to/site-packages --output /tmp/existing.json
python scripts/benchmark.py --installed /path/to/site-packages --output /tmp/timing.json --repeats 50
```

Hooks needed outside the owned files:

1. A public seeded equilibrium solve accepting a previous State's temperature,
   species amounts and active phases. `Solver::solve` always starts afresh;
   rocket now ports NASA's pressure Newton iteration and carries phase seeds,
   but cannot reuse full station guesses without this hook. tc_est currently
   uses public TP solves and equilibrium heat capacity, with verified enthalpy
   closure, to honor the estimate without private-core access.
2. Named Li(cr) reference data at 298.15 K, preserving its NASA formula, mass
   and zero reference enthalpy without exposing an invalid fit.
3. Equilibrium SP convergence/phase handling for example 13 at the exact state
   above. The failure is reproducible outside the rocket solver.

No equilibrium, reactant, package-version, CI or root README changes are
included. Transport at rocket stations was added later; see the transport row of
`docs/implementation-matrix.md`. NOTICE is unchanged. The rocket module and provenance file name
the Apache-2.0 NASA derivation.
