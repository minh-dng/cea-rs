# cea-rs

Experimental native Rust thermochemistry and chemical equilibrium, with Python
bindings. Canonical source: [minh-dng/cea-rs](https://github.com/minh-dng/cea-rs).
Install distribution `cea-keplr-rs`; import `cea_rs`. Version 0.1.0a3 is an
alpha, not a full NASA CEA replacement. Native `cea==3.3.4` remains the default
in the released Hybrid and keplr-solver consumers.

```python
from cea_rs import Solver

solver = Solver(["H2", "O2"])
state = solver.solve("TP", 3000.0, 1e6, [1.0, 8.0])  # K, Pa, mass weights
rocket = solver.rocket([1.0, 6.0], [298.15, 298.15], 3e6, 10.0,
                       flow="frozen_at_throat")
print(state.temperature_k, rocket.cstar_m_s, rocket.vacuum_velocity_m_s)
```

Implements TP/HP/SP/TV/UV/SV, experimental condensed/ion equilibrium,
IAC/FAC equilibrium or frozen rockets, incident/reflected shocks, equilibrium
CJ detonation, transport, and NASA input decks. Numerical work runs in
Rust, with no Fortran or native CEA runtime dependency. Python requires >=3.12.
The existing alpha packages remain in the original GitLab registry; this source
migration does not republish them. See [GitHub CI and releases](docs/github.md)
for the wheel targets and future release process.

[Implementation matrix](docs/implementation-matrix.md) gives exact units,
unsupported modes and failure behavior. [Validation](docs/validation.md)
separates the initial native comparisons, packaging tests and performance evidence.
[0.1.0a2 performance](docs/performance-0.1.0a2/README.md) records the reviewed
nozzle-root optimization; [historical CI validation](docs/ci-validation.md)
records local Linux checks and quota-blocked GitLab execution.
[GitHub Actions](https://github.com/minh-dng/cea-rs/actions) runs the current CI.
[Provenance](docs/provenance.md) identifies NASA and prior Rust sources/licenses.

Version 0.1.0a2 was **slower than native CEA** in its measured workloads.
Version 0.1.0a3 has not been benchmarked.
Full CEA parity and drop-in consumer compatibility are unproven. In particular,
N2O at Hybrid's 298.15 K input is rejected because its NASA fit starts at 300 K.

## Nozzle phase initialization

Equilibrium nozzle solves reuse the chamber's active condensed phases as
initial guesses. Phase removal, iteration budgets, equations, fitted-data and
closure checks remain unchanged. This fixes a carbon-bearing nozzle
initialization failure; it does not repair every standalone unseeded SP solve
or confer numerical acceptance. Consumers must update their exact version
guards and pins before loading alpha.3.

## Released consumers

- [Hybrid 3.3.0](https://gitlab.com/keplr-space/hybrid/-/blob/v3.3.0/docs/cea-rs-consumers.md)
  provides the `hybrid[cea-rs]` extra, pinned to `cea-keplr-rs==0.1.0a2`. Backend
  selection is explicit (`backend="cea_rs"`); the legacy 298.15 K N2O gas input
  also requires the documented `reference_298.15` policy, or an in-range inlet.
- [keplr-solver v0.1.8](https://gitlab.com/keplr-space/keplr-solver/-/blob/v0.1.8/docs/cea-rs-consumers.md)
  pins Hybrid 3.3.0 and exposes the same optional Rust package through
  `keplr-solver[cea-rs]`. Hybrid runs select `thermochemistry_backend="cea_rs"`;
  liquid requests select `requested_source="live_cea_rs"` under the documented
  policy constraints.

Installing either extra does not select Rust. These links describe released
library support, not Keplr application adoption or production deployment. The
consumer repositories and original packages still require GitLab access.
Scientific acceptance, full trajectories and NASA API parity remain unproven;
the consumer acceptance boundary remains `needs_human`.

## NASA input decks

```python
from pathlib import Path
from cea_rs import run_deck

cases = run_deck(Path("example11.inp").read_text())
rocket = cases[0].result
print(rocket.cstar_m_s, rocket.exit.state.temperature_k)
```

`run_deck(text)` returns a list of typed `DeckResult` objects. Each has
`case_name`, `case_index`, `kind`, `weights`, `result`, `transport`, and optional
`reactant_density_kg_m3`. The payload is an existing `State`, `RocketResult`,
`ShockResult`, or `DetonationResult`; all physical values use SI. Errors stop
the call without returning a partial list. Rust's `deck::parse_deck`,
`Deck::cases` and `Case::run` allow independent case diagnostics.

The parser covers NASA's sample `problem`, `reactants`, `only`/`omit`/`insert`,
`output` and `end` datasets, mixed units, custom formulas/enthalpies, mixture
ratios and station schedules. Unknown keywords or units fail by name.
Formatting/plot preferences are retained as metadata. Frozen detonation is
unsupported. No NASA report formatter, CLI or Python compatibility layer is
included. Existing solver numerics are unchanged.

All 55 cases in NASA's 14 main-interface decks solve. The table gives strict
`1e-4` relative comparison results against native 3.3.4 reports (example 6
uses NASA's complete historical report). Fractions >=1e-5 are compared;
zero-reference enthalpy has a 1e-3 J/kg absolute allowance.

| Example deck | Status | Worst relative difference |
|---|---|---:|
| 1: TP | Pass | 1.61392e-9 |
| 2: TV + transport | Pass | 6.74146e-8 |
| 3: HP | Pass | 4.40289e-10 |
| 4: UV | Pass | 1.39425e-10 |
| 5: custom binder | Pass | 1.48335e-5 |
| 6: detonation + transport | Pass | 4.11214e-5 |
| 7: equilibrium/frozen shock | Pass | 8.81075e-6 |
| 8: IAC schedules | Pass | 6.96275e-5 |
| 9: FAC contraction | Pass | 6.85920e-5 |
| 10: FAC mass flux | Pass | 6.68535e-5 |
| 11: ions + transport | Pass | 3.20473e-5 |
| 12: equilibrium + frozen | Fail: two trace fractions at area 75; nozzle stopping difference | 1.73452e-4 |
| 13: BeO insertion | Pass | 2.12241e-5 |
| 14: water condensation | Pass | 2.04574e-9 |
| Standalone sample 1 | Pass | 1.61392e-9 |
| Combined RP-1311 sample | 56/57 solve; known 1000 m/s shock failure and example 12 trace differences | 1.73455e-4 among successful cases |
| Frozen detonation | Unsupported; explicit error | — |

The unchanged historical reports are also tested. Examples 1, 7–11 contain
small historical differences above `1e-4`; example 12 covers only its frozen
branch. The suite records bounded exceptions for named cells and fails on
new differences. [Deck validation](docs/input-decks.md) gives both reference
tables, exact numbers, known NASA reporting/stopping limits and overlap with
existing hand-entered tests. [Fixture provenance](tests/data/nasa/README.md)
records NASA's notice, hashes and regeneration steps.

## Build

Use Rust 1.99.1, Python >=3.12 and maturin 1.15.0:

```sh
cargo +1.99.1 test --locked
python -m pip install maturin==1.15.0
maturin build --locked --release --out dist
```

The existing `mise.toml` is the maintainer's T7 checkout configuration and keeps
caches, temporary files and build outputs on that drive. Set `CARGO_HOME` after
mise resolves Rust. Use copy-mode Python installs on exFAT; its AppleDouble
sidecars can confuse tools that scan script directories.

This independent Keplr project is not endorsed by NASA. Apache-2.0; original
NASA and Bit Mango notices and license texts are retained. No consumer imports
or production defaults are changed by installing this separate package.
