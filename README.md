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
print(state.converged, state.temperature_k, rocket.complete, rocket.cstar_m_s)
```

Implements TP/HP/SP/TV/UV/SV, experimental condensed/ion equilibrium,
IAC/FAC equilibrium or frozen rockets, incident/reflected shocks, equilibrium
CJ detonation, transport, NASA input decks and a supported subset of NASA's
Python API in `cea_rs.nasa`. Numerical work runs in
Rust, with no Fortran or native CEA runtime dependency. Python requires >=3.12.
The existing alpha packages remain in the original GitLab registry; this source
migration does not republish them. See [GitHub CI and releases](docs/github.md)
for the wheel targets and future release process.

[Implementation matrix](docs/implementation-matrix.md) gives exact units,
unsupported modes and failure behavior. [Validation](docs/verbatim-parity.md)
records native comparisons, retained-state tests and remaining gaps.
[0.1.0a2 performance](docs/performance-0.1.0a2/README.md) records historical
measurements; [historical CI validation](docs/ci-validation.md)
records local Linux checks and quota-blocked GitLab execution.
[GitHub Actions](https://github.com/minh-dng/cea-rs/actions) runs the current CI.
[Provenance](docs/provenance.md) identifies NASA and prior Rust sources/licenses.

Version 0.1.0a2 was **slower than native CEA** in its measured workloads.
Version 0.1.0a3 has not been benchmarked.
The target is verbatim NASA CEA 3.3.4 behavior, including its quirks, failures
and false positives. Check returned convergence/status fields. A reported
success is not a proof of conservation or phase stability. Unsupported NASA
features raise NotImplementedError by name; see the implementation matrix.
The port does not call Fortran or the native wheel at runtime.

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
unsupported. No NASA report formatter or CLI is included. The numerical target is NASA CEA 3.3.4.

The 14 main-interface decks and combined RP-1311 sample are checked against
the pinned native reports. The 1000 m/s shock retains its native error code and
values; example 12 keeps the native area stop. Native false positives, frozen
workspaces and transport behavior are also covered by stored wheel snapshots.
[Current validation and remaining gaps](docs/verbatim-parity.md).
