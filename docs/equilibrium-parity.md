# Equilibrium parity

The target is NASA CEA 3.3.4 at commit `2f79a647e85737742fb08a127da00d1ea1d069da`.
NASA's Newton, truncation, phase and singular-recovery behavior determines the
answer and status. Independent conservation or stability checks do not replace
that behavior. Numerical failure retains the native fields.

`tests/fixtures/equilibrium-native.txt` and `tests/solver.rs` cover all six
modes, custom reactants, gas extensions, ions and condensed products.
`tests/test_nasa_api.py` checks mutable solution reuse and direct wheel output.
Condensed failure and false-positive snapshots are in
`tests/data/nasa/retained-3.3.4.json`.

The source was fetched into scratch space and built there. For the 300 K
Ar/NaCN case, its iteration matrices were compared with Rust to identify the
reference build's fused arithmetic and native species order. The production
solver does not invoke Fortran or import the native wheel.

[Condensed rules](condensed-phase-insertion.md), [current validation](verbatim-parity.md)
and [remaining gaps](implementation-matrix.md) define the current scope.
The earlier broad sweeps and performance JSON files remain historical evidence.
