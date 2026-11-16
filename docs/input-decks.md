# NASA input-deck support and validation

`cea_rs::deck::parse_deck(text)` reads one or more `end`-terminated problems
into `Vec<Deck>`. `Deck::cases()` expands schedules into typed `CaseInput`
values; `Case::run()` returns a `DeckResult`. This lets diagnostic tools
continue after a failed case. `cea_rs::deck::run_deck(text)` runs all cases
and returns `Result<Vec<DeckResult>>`, stopping at the first error without
returning a partial list.

Python exposes `cea_rs.run_deck(text) -> list[cea_rs.DeckResult]`:

```python
from pathlib import Path
from cea_rs import run_deck

for case in run_deck(Path("example11.inp").read_text()):
    print(case.case_name, case.case_index, case.kind)
    rocket = case.result  # RocketResult for a rocket deck
    print(rocket.cstar_m_s, rocket.exit.state.temperature_k)
```

`DeckResult.result` is `State`, `RocketResult`, `ShockResult` or
`DetonationResult`. `weights` follows input reactant order; `case_index` is
one-based within the named deck. `transport` contains the equilibrium/burned
state, all rocket stations, or shock stations 1, 2 and optionally 5, in that
order. Rocket `Station.transport` is also populated. Inlet densities, when
all are supplied, produce `reactant_density_kg_m3` by mass-weighted harmonic
averaging. Rust additionally retains the typed case inputs on each result.

## Supported datasets

| Dataset | Supported inputs |
|---|---|
| `problem` / `prob` | TP, HP, SP, TV, UV, SV; rocket/IAC/FAC; incident/reflected shock; equilibrium CJ detonation; case names, ions, equilibrium/frozen flags |
| State schedules | `p`, `t`, `rho`, `v`, `h/r`, `u/r`, `s/r`; all combinations for equilibrium; temperature × pressure for detonation |
| Mixture schedules | `o/f`, `f/o`, `f/a`, `%fuel`, mass `phi`, chemical `r` (including `r,eq.ratio`, `phi,eq.ratio`, `r,e`) |
| Rockets | `pi/p`, `pc/p`, `subar`/`sub,ae/at`, `supar`/`supersonic`/`sup-ae/at`, `acat`, `ma`/`mdot`, `nfz`/`nfrz`, assigned chamber `t` or `h/r`, `tcest` |
| Shocks | `u1` (m/s) or `mach1`; incident/reflected equilibrium/frozen combinations. Pressure/temperature schedules are paired with speed, repeating their final value, as NASA does |
| `reactants` / `react` / `reac` | `name`, `fuel`, `oxid`; mass or mole amounts, inlet temperature, density; exploded formulas, explicit molar enthalpy, optional molecular weight |
| `only`, `omit`, `insert` | Exact NASA product names and condensed seeds, passed to the public solver. Unknown omit names retain NASA's no-op behavior |
| `output` / `outp` | Transport and trace convergence; `massf`, `siunits`, `calories`, `long`, `short`, `debug`, and the samples' `plot` columns as retained report preferences |
| `end` | Required; datasets may appear in any order; several decks may share a file; `#` and `!` comments |

Keywords and units are case-insensitive; species names keep case and commas
(`C8H18(L),n-octa`). Supported problem keywords follow NASA's ordered
four-, three-, two- and one-letter prefix matching, rather than a list of
abbreviations: `rock`, `rocketry` and `rocket` all match `ro`, `shk` matches
`sh`, and `deton` matches `det`. Dataset names match their first four letters
(`end` uses three); reactant roles match two. Output keywords use NASA's
five-, four- and three-letter prefixes. Exploded formula elements without
a following numeric coefficient default to 1, including before the next
reactant parameter. These rules come from pinned `source/input.f90`,
lines 331–401, 497–530, 721–750, 807–829 and 987–1002.

The pinned parser accepts `tcesc`, `tcesr` and `tcesf`, but reads their
numeric value directly as the kelvin estimate (`input.f90:335`,
`main.f90:454`). Its unit-aware `parse_embedded` branch at `input.f90:872–877`
is bypassed for this parameter. The port follows that behavior:
`tcesc=3000` retains an estimate of 3000 K. Explicit delimited estimate units
are rejected rather than applying a conversion NASA never performs.

Numeric lists accept Fortran D exponents. Unknown keywords,
unknown units, missing values and unused physical constraints return named
errors. Mass and mole amounts, or `name` and `fuel`/`oxid`, cannot be mixed.
Each fuel/oxidant stream is normalized independently when a ratio is supplied.
Custom reference reactants use `Reactant::custom`; HP/rocket problems require
their enthalpy unless a chamber state is assigned. Formula-only reactants are
valid in temperature problems, where inlet enthalpy is not needed.

All parsed physical values and returned properties use the existing SI APIs.
Pressure constants match pinned NASA `units.f90`: 101325 Pa/atm,
6894.75729 Pa/psia and **133.3 Pa/mmHg**. Temperature conversions use
273.15 K and 1.8; 1 cal = 4.184 J. `h/r`, `u/r`, `s/r` multiply by
8314.51 on this port's kg basis. Specific-volume `v` schedules require an explicit unit,
as in NASA's `main.f90:1938–1945`; bare `v=...` is an error in TV, UV and SV.
Input volume in cc/g multiplies by 1e-3;
density in g/cc multiplies by 1000. Custom formula masses use NASA atomic
weights from the bundled neutral monatomic records. Report preferences never
change the SI result units or generate report/plot files.

Frozen detonation returns `Error::Unsupported` / `NotImplementedError`.
Frozen equilibrium state problems are also unsupported. Application solvers
retain their existing validation and reactant restrictions; for example,
custom reference enthalpies do not supply heat capacities for shocks.
Frozen rocket fit-limit stops retain `complete` and `stop_reason` from the
public API. The separate [NASA Python adapter](nasa-python-api.md) exposes
NASA calling conventions; there is no optional CLI.
No solver numerical code was changed.

## Reference policy

The fixtures include every problem `.inp` and example report found in the
pinned NASA repository and `cea==3.3.4` source distribution: 14 main-interface
decks/reports, the standalone sample of example 1, and the combined RP-1311
deck. Database `.inp` files are not problem decks. Originals are unedited;
NASA's notice and per-file SHA-256 hashes accompany the fixtures. Modern
reports were generated by an unmodified, independently rebuilt pinned NASA
executable; these reproduce the cached 3.3.4 executable byte-for-byte.
[Fixture provenance and regeneration](../tests/data/nasa/README.md).

The small report reader lives only in `tests/nasa_decks.rs`. It compares T,
P, density, H, U, S, M, Cp, gamma, sound speed, and mole/mass fractions at
or above 1e-5. Rockets add area, Mach, c-star, Cf, Isp and vacuum Isp;
shocks add initial/downstream velocities and jump ratios; detonations add
unburned properties, wave speed, Mach and jump ratios. Every available
transport Cp/viscosity/conductivity/Prandtl row is compared. Imperial c-star
and impulse rows are converted to SI. All nonzero values use **1e-4 relative
tolerance**. Zero Cp uses 1e-4 absolute tolerance; NASA's zero reference
unburned enthalpy uses 1e-3 J/kg absolute tolerance (the fit gives
-0.0005062602429 J/kg). Zero-reference enthalpy has no relative difference
and is excluded from the reported worst relative error.

Tests check every successful case and every reported rocket or shock station,
including repeated report-page prefixes. The virtual FAC infinity station
is absent from NASA's printed tables. Known differences are allowed only
for explicitly named cells with bounded errors; every other cell keeps
the strict tolerance. A new discrepancy fails the test. A green regression
suite therefore does **not** mean every historical/reference cell passes.

## Per-deck results

“Primary” uses the modern 3.3.4 report, except example 6, which uses the
complete historical long report because of the modern printer defect below.
All 55 main-interface cases solve successfully. Differences are relative.

| Deck | Cases | Primary status / worst difference | Shipped historical report |
|---|---:|---|---|
| 1, TP H2/Air | 12 | Pass, 1.61392e-9 | Fail, 1.10248e-4 in two first-case carbon fractions |
| 2, TV + transport | 3 | Pass, 6.74146e-8 | Pass, 8.18672e-6 |
| 3, HP liquid fuels | 3 | Pass, 4.40289e-10 | Pass, 2.23315e-8 |
| 4, UV | 1 | Pass, 1.39425e-10 | Pass, 3.40059e-8 |
| 5, custom binder + condensates | 5 | Pass, 1.48335e-5 | Pass, 1.65275e-5 |
| 6, detonation + transport | 4 | Pass, 4.11214e-5 | Same reference; modern report unsuitable |
| 7, equilibrium/frozen incident shock | 12 | Pass, 8.81075e-6 | Fail, 1.68507e-4; old mmHg conversion |
| 8, IAC schedules | 1 | Pass, 6.96275e-5 | Fail, 1.14411e-4; subsonic Mach/Cf/Isp |
| 9, FAC contraction ratio | 1 | Pass, 6.85920e-5 | Fail, 1.14472e-4; combustor Mach/Cf/Isp |
| 10, FAC mass flux | 1 | Pass, 6.68535e-5 | Fail, 1.39245e-4; combustor area/Mach/Cf/Isp |
| 11, ions + station transport | 1 | Pass, 3.20473e-5 | Fail, 1.00125e-4; subsonic Mach/Cf/Isp |
| 12, equilibrium + throat-frozen | 2 | **Fail**, 1.73452e-4; two equilibrium trace mass fractions | Frozen branch passes, 7.11817e-6; equilibrium branch absent |
| 13, inserted BeO phase | 1 | Pass, 2.12241e-5 | Pass, 1.93431e-5 |
| 14, water condensation | 8 | Pass, 2.04574e-9 | Pass, 2.99496e-7 |
| `samples/example1.inp` | 12 | Pass, 1.61392e-9 | Same physical case as 1 |
| `samples/rp1311_examples.inp` | 57 | 56 solve; one known failed shock. Successful cases compared; example 12 has the same two trace mismatches, worst 1.73455e-4 | Combined deck has no original standalone report |
| `shock-reflected-transport.inp` | 4 | **Fail**, 2.06110e-2; two frozen reflected header Cp cells; all other state, composition, ratio and transport cells pass | Additional regression deck; modern report only |
| Frozen detonation (validation input) | — | Unsupported, explicit error | Unsupported by NASA 3.3.4 |

The combined deck's valid shock cases pass at 9.19966e-6. Its seventh and
twelfth blocks have their own modern reports because the pressures/speeds
and HNO2 product pool differ from the main-interface decks. Only `output long`
was added when obtaining those reports; no physical input was changed.
The other combined blocks use the equivalent per-example physical models.

## Solver and reference discrepancies

No new port equilibrium defect was found. The strict example 12 nozzle
comparison remains a failure; solver stopping and reference differences are
recorded rather than repaired on this branch.

1. **Example 12, equilibrium exit at requested area 75.** Rust's visible
   station 7 (eighth station) has area 75.00000000085697, T=1439.5374600082164 K,
   P=7587.699821861046 Pa. NASA reports area 74.99733853978925,
   T=1439.549820428637 K, P=7588.049271266402 Pa. Rust/NASA mass fractions are
   CO=4.900777497427966e-5 / 4.901627693862203e-5 (1.73452e-4 relative) and
   OH=3.417414926408429e-5 / 3.417777936134580e-5 (1.06212e-4).
   At NASA's T/P, public Rust TP gives CO=4.901627662510234e-5 and
   OH=3.417777946875311e-5, agreeing within 6.40e-9 relative. The test checks
   this independently and requires Rust's area ratio to be within 1e-8
   of 75 for both main-interface and combined inputs. NASA's supersonic
   equilibrium area tolerance is 4e-5 (`source/rocket.f90:870,920`).
   NASA's residual area error is -3.54861e-5 relative;
   the small state difference amplifies in trace composition. The combined
   deck's additional HNO2 changes CO to 4.900914902465053e-5 (Rust) versus
   4.901765137682429e-5 (NASA), a 1.73455e-4 difference at the same exit.

2. **Combined example 7, equilibrium shock at 1000 m/s.** Inputs are
   H2/O2/Ar moles 0.05/0.05/0.9, 300 K, 10 mmHg = 1333 Pa. Rust returns
   `not_converged: case "7" #1: Shock exceeded 60 iterations`.
   NASA also prints `WARNING: No convergence` and `ANSWERS NOT RELIABLE`,
   but retains P2=8870.053226803845 Pa, T2=1371.901376071067 K and
   u2=703.5613575871638 m/s. Its printed T2/T1=4.656334035059522 disagrees
   with the retained T2/300=4.573004586903557. That column is not an oracle;
   all 13 other shock cases are compared at their original report columns.
   The test asserts the exact Rust error and keeps a failed-case slot,
   without removing any reference values. This overlaps the existing
   `native_nonconvergence_never_returns_a_retained_incident_or_reflected_state`
   regression. Nothing stale is returned by this port.

3. **Modern example 6 printer.** The section labeled `UNBURNED GAS` reports
   burned T=3674.280831065507 K for a 298.15 K inlet. The `BURNED GAS` section
   then repeats rounded values, insufficient for some strict comparisons.
   The original long report contains correct unburned and burned tables and
   is used directly, including transport and jump ratios.

4. **Historical references.** Example 1's first CO fraction is
   1.707346998959570e-4 versus Rust 1.707158766768762e-4; CO2 has the same
   1.10248e-4 relative drift. The modern report agrees with Rust. Example 7
   uses 1333.223684210530 Pa for 10 mmHg, whereas pinned `units.f90` uses
   1333 Pa; pressure, density and a few H fractions consequently exceed
   1e-4. Examples 8–11 differ only in named subsonic/combustor performance
   cells. Their worst Mach pairs (Rust / historical NASA) are
   0.4131280218181437 / 0.4131752938039900 (8),
   0.4131488008983451 / 0.4131961002646020 (9),
   0.4128501072729505 / 0.4129076025319600 (10), and
   0.05941049944432838 / 0.05941644852028900 (11).
   Example 10's combustor area is 1.581084992772664 / 1.580904539442940.
   All pass against the modern executable. Example 12's original report
   requests frozen only; it does not validate the new equilibrium branch.

5. **Frozen reflected header Cp.** The additional transport deck uses
   example 7's H2/O2/Ar mixture, 300 K, 10 mmHg and 1400 m/s, with all four
   incident/reflected equilibrium/frozen combinations. NASA prints stations
   1/2 for the first and third cases, and station 5 for every case. Its
   equilibrium-reflection sections do not repeat stations 1/2
   (`source/main.f90:1046–1073`). Transport is printed for stations 2 and 5,
   with only frozen transport rows for frozen shocks
   (`main.f90:1220–1222,1358–1359,1474–1516`); no initial transport rows
   are printed. All 187 eligible cells are checked, with only two bounded
   exceptions: case 1 station 5 Cp is 602.9576961140284 J/(kg K) in Rust
   versus 593.2530916177532 in NASA (1.63583e-2 relative, bound 1.64e-2);
   case 3 is 603.5207357881964 versus 591.3327704887254
   (2.06110e-2, bound 2.07e-2). These are physical Cp at T5 versus NASA's
   carried incident Cp header basis, explicitly assigned in
   `source/shock.f90:374–377`. NASA's separate frozen transport Cp uses T5
   and agrees with the port at 1e-4. Every other newly compared cell passes
   1e-4; no solver numerics or physical Cp fields were changed.

## Overlap with existing hand-entered tests

These tests remain intact; the deck suite replaces no checks:

- `native_equilibrium_domain_regressions`: examples 1–5 and 14 from
  `tests/fixtures/equilibrium-native.txt`.
- `matches_native_cj_states_and_conserves_fluxes`: example 6.
- `native_shock_goldens_and_conservation`: example 7's incident modes, with
  the explicitly documented historical/modern mmHg difference.
- `rocket_stations_match_independent_native_goldens`: examples 8–12; its
  older example 11 fixture uses a Li reference, whereas the deck uses Li(cr).
- `rp1311_example_13_continues_through_beo_coexistence`: example 13 stations.
- `nasa_transport_goldens`: example 2 TV states, example 6 burned states and
  example 11 chamber TP; `rp1311_example_11_station_transport_matches_native`
  overlaps the full example 11 station/transport comparison.

These older tests also check conservation, failure behavior and extra native
cases that the example reports do not cover. The input-deck work does not
replace those purposes.

## Build and wheel checks

Validated on macOS arm64, Rust 1.99.1, Python 3.15.1, 2026-11-17, starting
from `300328e`. Build output, temporary files and uv cache used local temporary
directories; `mise.toml` was not edited.

| Check | Result |
|---|---|
| Baseline `cargo test --release` | 62 passed, 1 ignored |
| Final `cargo test --release` | 68 passed, 1 ignored; existing ignored diagnostic retained |
| `cargo clippy --release --all-targets -- -D warnings` | Exit 0 |
| `cargo clippy --release --all-targets --features python -- -D warnings` | Exit 0 |
| `cargo fmt --check` | Exit 0 |
| `uvx maturin build --release` | Built `cea_keplr_rs-0.1.0a3-cp312-abi3-macosx_11_0_arm64.whl` |
| Installed wheel, `tests/test_wheel.py` | 6 passed, including all application result types and deck errors |
| NASA original bytes and `SHA256SUMS` | All 47 source/reference files verified |

Commands (the wheel test uses mise's selected Python directly):

```sh
S=$(mktemp -d)
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo test --release
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets --features python -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo fmt --check
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" UV_CACHE_DIR="$S/uv-cache" \
  uvx maturin build --release --out "$S/wheels"
mise exec -- env TMPDIR="$S" UV_CACHE_DIR="$S/uv-cache" \
  uv pip install --python "$(mise which python)" --target "$S/installed" "$S"/wheels/*.whl
mise exec -- python tests/test_wheel.py --installed "$S/installed"
```

On this machine `mise exec -- env ... python` selects global Python 3.10,
which cannot load the Python >=3.12 wheel; `mise exec -- python` selects the
configured 3.15.1. The initial import failure was reproduced and resolved by
using that interpreter; no packaging or solver changes were needed.
