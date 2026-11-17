# cea-rs

Chemical equilibrium and rocket performance in Rust, with Python bindings.
It copies the behaviour of [NASA CEA](https://github.com/nasa/cea) 3.3.4,
including its quirks and failure codes, and needs no Fortran at runtime.
Version 0.1.0a4 is an alpha.

## Install

Python 3.12 or later. Building from source needs a Rust toolchain
([rustup](https://rustup.rs)).

```sh
uv add git+https://github.com/minh-dng/cea-rs
# or
pip install git+https://github.com/minh-dng/cea-rs
```

Prebuilt wheels for Linux and macOS are attached to
[GitHub releases](https://github.com/minh-dng/cea-rs/releases) when published.

## Use it in place of NASA CEA

`cea_rs.nasa` has the same classes, arguments and units as NASA's `cea`
package, so most scripts only change the import:

```python
import numpy as np
import cea_rs.nasa as cea  # was: import cea

reac = cea.Mixture(["H2", "O2"])
prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
solver = cea.EqSolver(prod, reactants=reac)
soln = cea.EqSolution(solver)
solver.solve(soln, cea.TP, 3000.0, 10.0, reac.moles_to_weights(np.array([2.0, 1.0])))
print(soln.converged, soln.T, soln.mole_fractions["H2O"])
```

As in NASA CEA, a failed solve does not raise. It sets `converged=False` and
`last_error` and gives a `RuntimeWarning`, so check `converged`. Features not
yet ported raise `NotImplementedError` by name.

## SI API

```python
from cea_rs import Solver

solver = Solver(["H2", "O2"])
state = solver.solve("TP", 3000.0, 1e6, [1.0, 8.0])  # K, Pa, mass weights
rocket = solver.rocket([1.0, 6.0], [298.15, 298.15], 3e6, 10.0,
                       flow="frozen_at_throat")
print(state.converged, state.temperature_k, rocket.complete, rocket.cstar_m_s)
```

It covers TP/HP/SP/TV/UV/SV equilibrium, condensed phases and ions, IAC/FAC
rockets, incident and reflected shocks, CJ detonation, transport properties,
and NASA input decks (see below).

## Documentation

- **NASA CEA documentation:** <https://nasa.github.io/cea/>. Use it for the
  `cea_rs.nasa` API, input decks and the underlying theory.
- [NASA Python API inventory](docs/nasa-python-api.md): every NASA class and
  member, its units, and whether it is supported here.
- [Implementation matrix](docs/implementation-matrix.md): what this port
  supports, with units and failure behaviour.
- [Parity with NASA](docs/verbatim-parity.md): native comparisons, the
  largest differences and the remaining gaps.
- [Future improvements](docs/future-improvements.md): corrections to NASA
  behaviour that were removed and may return as options.
- [Input decks](docs/input-decks.md), [CI and releases](docs/github.md),
  [provenance and licences](docs/provenance.md).

Version 0.1.0a2 was slower than native CEA in its measured workloads
([measurements](docs/performance-0.1.0a2/README.md)); 0.1.0a4 has not been
benchmarked.

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
