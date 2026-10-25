# Validation evidence for 0.1.0a1

Scientific comparisons use independently installed native NASA `cea==3.3.4`,
matching the existing Hybrid seam, under its Python 3.15.0 interpreter. `docs/native-comparison.json` records the
loaded oracle shared-object/data hashes, Rust data hash, environment, every
native/Rust result and difference. Each of nine fixture groups runs in a
separate subprocess; native aborts/crashes/timeouts fail that group and the
merged report. Native output never enters the Rust computation.

The fixed set has 146 comparison rows: 16 equilibrium states in six modes,
six composition checks, four condensed/ion checks, and 120 rocket station/
performance checks from 30 complete rockets. Cases include H2/O2,
CH4/O2/N2, graphite from CO, a declared H/H2/H+/e- ion pool, LOX/RP-1,
N2O/paraffin, custom N2O/ABS, and NTO/Aerozine50. Rockets use 1 MPa/AeAt=3.68
and 5 MPa/AeAt=20, in equilibrium, chamber-frozen and throat-frozen flow.

The difference is `abs(Rust-native)/max(abs(native),1 SI unit)`. Assertion
limits: 3e-4 for equilibrium properties, 3e-5 absolute composition difference,
and 1e-4 for rocket properties/velocities/c-star. Observed maximum rocket
station difference is approximately 5.84e-5 (0.00584%); this is a finite-case
comparison, not a general accuracy guarantee. Frozen equilibrium Cp is not
compared: the Rust API correctly marks it unavailable.

The earliest manually assembled ion fixture omitted ions=True on the native
reactant mixture and native CEA aborted. That fixture was corrected; the
released comparison runner process-isolates the fixtures. This native crash
is not counted as a solver pass. Canonical N2O at 298.15 K remains a tested
DomainError, not a passing compatibility case.

Rust checks: 18 inherited asserting unit tests (one upstream diagnostic sweep
ignored because it does not assert success) plus five new integration tests.
The integration tests check six-mode round trips, rejection of invalid/
impossible problems, graphite and ions, and sonic/mass/entropy/energy closure
for all three nozzle modes. The vendored thermodynamic crate also has tests.

Installed-wheel checks run with `python -I` from outside the source tree under
Python 3.12.14 and 3.15.1. They verify import provenance, invalid/domain/
unsupported failures, low-temperature TP/hot-rocket/frozen interleaving,
interleaved NTO/Aerozine50 calls, concurrent calls, and nozzle identities.
Repeated Rust results are bitwise identical in those checks. This is not a
Hybrid cold-flow/hot-fire simulation acceptance claim.

```sh
# Use an interpreter with the independent native cea==3.3.4 + numpy installed.
python scripts/compare_native.py --installed /external/fresh-wheel-site --output docs/native-comparison.json
python scripts/benchmark.py --installed /external/fresh-wheel-site --output docs/benchmark.json
# Run outside the repository, passing an absolute script path.
python -I /absolute/cea-rs/tests/test_wheel.py --installed /external/fresh-wheel-site
```

The release-mode benchmark records construction, first-call, and 50 repeated
warm solves after five discarded warmups, with a fresh native solution object
per call as in Hybrid. Native/Rust measurement order alternates. OS disk
caches are not flushed; first-call figures are process-local, not disk-cold
measurements. An additional benchmark launches five fresh Python processes
per backend for TP, equilibrium rockets and throat-frozen rockets, recording
import/init, construction and first solve separately. The OS filesystem cache
is not flushed; this is process-cold, not disk-cold. Validation stays enabled in Rust. `benchmark.json` records raw
median/min/max and warm ratios. This implementation is slower than native:
repeated independent equilibrium solves in bracketed nozzle roots dominate.
No production throughput or full simulation benchmark was run.

Artifact inspection verifies wheel/sdist license notices and rejects AppleDouble
`._*` entries. Maturin exclusions remove those filesystem sidecars, including
matches from license-file globs. The rebuilt extension's dynamic dependencies
are macOS system libraries; no libgfortran or CEA library is linked. Sdist builds
and fresh installs use separate external directories, not an editable install.
