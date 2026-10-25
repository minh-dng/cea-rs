# cea-rs

Experimental native Rust thermochemistry and chemical equilibrium, with Python
bindings. Install distribution `cea-keplr-rs`; import `cea_rs`. Version 0.1.0a1
is an alpha, not a full NASA CEA replacement. Native NASA CEA remains the Keplr
and Hybrid production default.

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
Rust, with no Fortran or native CEA runtime dependency. Python requires >=3.12;
the first binary wheel targets macOS arm64, with an sdist for other builds.

[Implementation matrix](docs/implementation-matrix.md) gives exact units,
unsupported modes and failure behavior. [Validation](docs/validation.md)
separates native comparisons, packaging tests and performance evidence.
[Provenance](docs/provenance.md) identifies NASA and prior Rust sources/licenses.

The current release is **slower than native CEA** in the measured workloads.
Full CEA parity and drop-in consumer compatibility are unproven. In particular,
N2O at Hybrid's 298.15 K input is rejected because its NASA fit starts at 300 K.

Build using the existing mise toolchain. This checkout deliberately keeps
Cargo/uv caches, temporary files and outputs on the T7 external drive:

```sh
mise trust
mise exec -- env CARGO_HOME="/Volumes/T7 Shield/keplr-work/cache/cargo" cargo test --locked
mise exec -- env CARGO_HOME="/Volumes/T7 Shield/keplr-work/cache/cargo" maturin build --locked --release --out dist
```

Set Cargo home after mise resolves the toolchain: overriding `CARGO_HOME` before
mise runs can make its rustup detection repeatedly attempt installation. Use
copy-mode Python installs on exFAT; its AppleDouble sidecars can confuse tools
that scan script directories. Ruff checks ignore EXE002 only because exFAT
marks regular source files executable.

This independent Keplr project is not endorsed by NASA. Apache-2.0; original
NASA and Bit Mango notices and license texts are retained. No consumer imports
or production defaults are changed by installing this separate package.
