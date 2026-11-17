# NASA CEA 3.3.4 parity

The target is the behavior of NASA CEA at commit
[`2f79a647e85737742fb08a127da00d1ea1d069da`](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da).
The independent numerical reference is the PyPI `cea==3.3.4` wheel.
A better physical answer, a stricter acceptance check or an extra recovery
solve is outside this target. The port keeps NASA's quirks and reported failures.

The pinned Fortran source was fetched and built outside the repository.
Detailed Fortran/Rust matrix traces identified the native species order and
fused arithmetic needed for the 300 K NaCN case. Rust does not invoke Fortran
or import the reference wheel when solving. Bundled NASA data are unchanged.

## Behavior covered

The condensed solver follows native insertion, removal, phase numbering,
coexistence, singular recovery, pivot limits, SMNOL and iteration limits.
Cold and fixed-temperature recovery solves have been removed.

Rust states expose `converged` and `last_error`. The NASA Python adapter warns
on numerical failure, retains the fields NASA retains, and returns normally.
A reused failed solution keeps the previous report and transport fields while
exposing the new raw temperature and species amounts. Positive fixed
temperatures outside 160–22000 K return NASA's post-processed code 8 state.

The tests preserve these native results:

- Ar/Hg HP and the water/NaCN transition failures, and NaCN success at 300 K.
- Br2 supersaturation without solid insertion.
- Air/Jet-A at 20 MPa reporting success despite conservation defects.
- AP/Al at 300 K reporting success with NaN thermodynamic fields and zero transport on a fresh solution; reused solutions retain prior transport.
- The 1000 m/s RP-1311 shock retaining code 9; its reflected calculation still continues.
- Example 12's stop near area 74.99733854 for requested area 75.
- FAC early freezing, including partial results and NaN workspaces.
- Native frozen Cp/Cv fields and transport Cp replacement.

Shock and detonation use NASA's analytic equations. Rocket throat, area,
frozen and FAC iterations use native stopping tests. Rocket performance
recomputes chamber enthalpy from stored composition, separately from the flow
iteration. Transport uses the gas block, fills native missing-data estimates,
and uses the overwritten reaction-system right-hand side after singularity.

Every removed correction is listed in [future improvements](future-improvements.md).
There are no flags that restore those corrections.

## Validation on macOS arm64

Checks use Rust 1.99.1, CPython 3.12.14 and NumPy 2.5.5.

| Check | Result |
|---|---|
| Release Rust suite | 61 passed; one existing diagnostic ignored |
| Installed-wheel tests | 8 passed |
| NASA Python tests with the native wheel | 9 passed, including all flag and schedule subcases |
| NASA Python tests without the native package | 9 passed against independent snapshots |
| Broad native sweep | 146 comparison rows passed |
| Rocket schedule sweep | 24 cases passed, including matched nonfinite fields and statuses |
| Clippy | Default and Python-feature all-target checks passed with warnings denied |
| Formatting | `cargo fmt --check` passed |

The Python tests use 41 successful native snapshots and 13 retained-state
snapshots. Live comparisons run the same calls against the independent wheel.
The Rust deck tests cover all 14 main-interface decks and the combined RP-1311
sample. Modern printed reports have no per-cell tolerance exemptions.
Historical CEA2 reports retain explicitly bounded printing/version differences.

Across the unrounded native comparisons, the worst relative difference is
**7.554308465184085e-5 (0.0075543%)**. It is the chamber enthalpy in the H2/O2,
5 MPa, area-20, throat-frozen case: native `-0.0005354430331850154 J/kg`,
Rust `-0.0005354025841666333 J/kg`. The absolute difference is about
`4.04e-8 J/kg`. The next largest difference is about `6.52e-5` in a FAC
throat-frozen performance coefficient close to zero. The NASA Python suite's
worst difference is `3.712466258753133e-6`, example 11 equilibrium conductivity.
The historical 2004 printed reports have a larger maximum, `1.67777e-4`, from
example 7's old pressure conversion; they are not the 3.3.4 numerical oracle.

[Native sweep](native-comparison.json), [rocket sweep](rocket-parity/comparison.json)
and [Python results](nasa-python-api-results.json) contain the numerical evidence.
Relative differences use nonzero native values; zero fields use absolute bounds.
NaN locations and status codes must match. Some older sweep assertions use
`max(abs(native), 1 SI unit)` scaling; the maximum above was recomputed from
the recorded unscaled values. Trace species and finite-difference noise use
the explicit bounds in their tests.

This establishes parity for the compared inputs. It does not establish
bit-for-bit agreement for every chemistry or compiler.

## Remaining gaps

The [implementation matrix](implementation-matrix.md) and
[Python member inventory](nasa-python-api.md) list the supported surface.
Unported Python features retain their names and arguments and raise a named
`NotImplementedError`:

- Smooth truncation, custom binary databases and fitted-reactant overrides.
- Frozen-reactant rockets, frozen detonation and frozen rocket gas logs.
- Some ion, reduced-element, condensed UV/SV and absent-element sensitivities;
  derivatives of nonconverged solutions and verbose derivative diagnostics.
- MATLAB wrappers and diagnostic logging other than `LOG_NONE`.
- Unsolved flow report fields and uncomputed derivative workspaces. Native
  flow storage is uninitialized; reading it can expose garbage or crash.
- Single-product pools with incompatible element inventories. The reference
  wheel crashes on the H2/O2 1:1, H2O-only case, so this is a named gap.

EqSolution initial guesses are supported, including NASA's Python setter
quirk: the first solve restores the saved constructor seed. NASA's broken
shock amount getters still raise its missing-`_get_size` AttributeError.
There is no NASA CLI, report/plot generator, custom database loader or
C/Fortran/Excel ABI. Linux/Windows rounding and production consumers were not
validated in this pass.

## Reproduce the checks

Keep build files outside the repository. The drive in `mise.toml` need not be
mounted, and that file must not be changed.

```sh
S=$(mktemp -d /tmp/cea-parity.XXXXXX)
mkdir -p "$S/tmp"
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S/tmp" cargo test --release
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S/tmp" cargo clippy --all-targets -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S/tmp" cargo clippy --all-targets --features python -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S/tmp" cargo fmt --check
```

Create a separate Python environment, build a release wheel, and install it
outside the source package. Set `UV_CACHE_DIR="$S/uv-cache"` and
`TMPDIR="$S/tmp"` when using uv through mise.

```sh
uv venv --python 3.12 "$S/venv"
uv pip install --python "$S/venv/bin/python" maturin numpy cea==3.3.4
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S/tmp" \
  "$S/venv/bin/maturin" build --release --out "$S/wheels"
uv pip install --python "$S/venv/bin/python" --target "$S/installed" "$S"/wheels/*.whl
"$S/venv/bin/python" tests/test_wheel.py --installed "$S/installed"
"$S/venv/bin/python" tests/test_nasa_api.py --installed "$S/installed" --report "$S/python.json"
"$S/venv/bin/python" scripts/compare_native.py --installed "$S/installed" --output "$S/native.json"
PYTHONPATH="$S/installed" "$S/venv/bin/python" scripts/compare_rocket_parity.py --output "$S/rocket.json"
"$S/venv/bin/python" scripts/record_retained_states.py --installed "$S/installed" --output "$S/retained.json"
```

Repeat the Python tests in an environment containing NumPy and the candidate
wheel, without installing `cea`, to check the stored snapshots independently.
[Fixture instructions](../tests/data/nasa/README.md) explain the Fortran reports
and native snapshot generation. No native expected value comes from a Rust solve.
