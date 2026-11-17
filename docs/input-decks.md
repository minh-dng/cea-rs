# NASA input-deck support and validation

`cea_rs::deck::parse_deck(text)` reads one or more `end`-terminated problems
into `Vec<Deck>`. `Deck::cases()` expands schedules into typed `CaseInput`
values; `Case::run()` returns a `DeckResult`. Numerically failed cases retain their native state and status. `cea_rs::deck::run_deck(text)` runs all cases
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
keep the named unsupported-feature limits in the implementation matrix; for example,
custom reference enthalpies do not supply heat capacities for shocks.
Frozen rocket fit-limit stops retain `complete` and `stop_reason` from the
public API. The separate [NASA Python adapter](nasa-python-api.md) exposes
NASA calling conventions; there is no optional CLI.
Solver acceptance and stopping behavior follow NASA 3.3.4, including false positives.

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

Tests check every available case and every reported rocket or shock station,
including repeated report-page prefixes. The virtual FAC infinity station
is absent from NASA's printed tables. Known differences are allowed only
for explicitly named cells with bounded errors; every other cell keeps
the strict tolerance. A new discrepancy fails the test. A green regression
suite therefore does **not** mean every historical/reference cell passes.

## Native failures and reporting quirks

All main-interface samples and the combined RP-1311 cases return their native
states. The combined 1000 m/s shock carries code 9 and its retained fields.
Example 12 keeps the native area stop near 74.99733854 for requested area 75;
its CO/OH cells no longer have a separate tolerance. Reflected-frozen shock
header Cp retains NASA's incident-state basis, while transport Cp uses the
current gas mixture. These are distinct stored fields.

The new native 3.3.4 reports are the primary oracle. Original RP-1311 reports
remain historical fixtures: their old unit rounding and printed values have
explicitly bounded differences. They do not override the pinned source or wheel.
The modern example-6 printer mislabels the unburned table, so its complete
historical report and native Python snapshots provide complementary checks.

The parser supports the keywords listed above. It does not generate a NASA
report, plot, CLI or ABI. Unknown physical keywords fail by name. Frozen
detonation and other unported physics raise a named unsupported-feature error.
See [current checks](verbatim-parity.md) and [remaining gaps](implementation-matrix.md).
