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

Implements TP/HP/SP/TV/UV/SV, experimental condensed/ion equilibrium, and
infinite-area equilibrium or frozen rocket expansion. Numerical work runs in
Rust, with no Fortran or native CEA runtime dependency. Python requires >=3.12.
The existing alpha packages remain in the original GitLab registry; this source
migration does not republish them. See [GitHub CI and releases](docs/github.md)
for the wheel targets and future release process.

[Implementation matrix](docs/implementation-matrix.md) gives exact units,
unsupported modes and failure behavior. [Validation](docs/validation.md)
separates the initial native comparisons, packaging tests and performance evidence.
[Alpha2 performance](docs/performance-0.1.0a2/README.md) records the reviewed
nozzle-root optimization; [historical CI validation](docs/ci-validation.md)
records local Linux checks and quota-blocked GitLab execution.
[GitHub Actions](https://github.com/minh-dng/cea-rs/actions) runs the current CI.
[Provenance](docs/provenance.md) identifies NASA and prior Rust sources/licenses.

Alpha.2 was **slower than native CEA** in its measured workloads; alpha.3 has
no new performance claim.
Full CEA parity and drop-in consumer compatibility are unproven. In particular,
N2O at Hybrid's 298.15 K input is rejected because its NASA fit starts at 300 K.

## Alpha.3 nozzle continuation

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
