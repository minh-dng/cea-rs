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
Repeated-singularity retention requires a reduced element set and positive
temperature and gas moles. Condensed-set convergence exhaustion leaves the
earlier report alone. Retained states use NASA's post-processing without new
partials, and transport seeds are restored on the same exits as in NASA.

Rust states expose `converged` and `last_error`. The NASA Python adapter warns
on numerical nonconvergence, retains the fields NASA retains, and returns normally.
A Fortran abort raises `RuntimeError` after updating the solution. The abort
keeps NASA's `converged=True`, code 10 and working iterate; the next solve
continues from that iterate.
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
- Shock gas moles initialized as `1/sum(weights)`, including unnormalized weights.
- N2 ion polishing at 300 K retaining its electron correction when the charge
  sum underflows, then returning code 8 and `n=4311816.4723374555`.
- Zero-weight reactants retaining their elements in the product pool.
- Ar/water abort/reuse and Br2 transport-seed restoration across no-gas states.
- Failed equilibrium shocks retaining volume after density is cleared; failed
  frozen shocks keep their unfilled volume at zero.
- Mach-1.01 air with frozen reflection retaining NASA's singular incident
  conductivity and Prandtl number. The shock species amounts now use
  `(weight/molecular_weight)/sum(weights)`, in NASA's operation order.
- Failed equilibrium reflections and frozen temperature-cap exits with cleared
  report fields and zero shock ratios.
- Mercury incident shocks with no valid state at 100 K and a retained last
  valid state at 500 K, both at Mach 1.01 and 1 bar.
- The atomic transport basis for an equilibrium incident shock followed by a
  frozen reflection, including a tighter air-conductivity comparison.
- Mercury and hot H2/O2 detonations that continue with nonfinite internal
  iterates and return status 8 rather than raising an input-domain exception.

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
| Release Rust suite | 62 passed; one existing diagnostic ignored |
| Installed-wheel tests | 7 passed with native CEA; 7 passed without it |
| NASA Python tests with the native wheel | 10 passed; one explicit expected numerical failure |
| NASA Python tests without the native package | 10 passed; the same expected numerical failure against an independent snapshot |
| Broad native sweep | 146 comparison rows passed |
| Rocket schedule sweep | 24 cases passed, including matched nonfinite fields and statuses |
| Clippy | Default and Python-feature all-target checks passed with warnings denied |
| Formatting | `cargo fmt --check` passed |

The Python tests use 41 successful native snapshots and 31 retained-state
snapshots. Every retained-state snapshot now contains numerical fields.
The singular ion-shock comparison of P21, T21, station thermodynamics, flow
and transport is an explicit `expectedFailure`, named
`singular-ion-shock-roundoff` below. Its status and positive finite station
fields are checked separately by an ordinary passing test. Live comparisons
run the same calls against the independent wheel.
The Rust deck tests cover all 14 main-interface decks and the combined RP-1311
sample. Modern printed reports have no per-cell tolerance exemptions.
Historical CEA2 reports retain explicitly bounded printing/version differences.

Across the passing numerical comparisons, excluding the explicitly expected
`singular-ion-shock-roundoff` failure, the worst relative difference is
**7.554308465184085e-5 (0.0075543%)**. It is the chamber enthalpy in the H2/O2,
5 MPa, area-20 schedule: native `-0.0005354430331850154 J/kg`,
Rust `-0.0005354025841666333 J/kg`. The absolute difference is about
`4.04e-8 J/kg`. The next largest difference is about `6.52e-5` in a FAC
throat-frozen performance coefficient close to zero. The NASA Python suite's
worst difference is `3.712466258753133e-6`, example 11 equilibrium conductivity.
The historical 2004 printed reports have a larger maximum, `1.67777e-4`, from
example 7's old pressure conversion; they are not the 3.3.4 numerical oracle.

[Review reproductions](retained-state-review.json) record the seven additional
before/after comparisons, including the air transport inputs. No new numerical
gap or mixture-specific correction was needed. The supplied `reproduce.py`
also passes a full-field comparison at 1e-9 relative tolerance. The final
broad and rocket numerical records are unchanged from the linked sweeps.

[Native sweep](native-comparison.json), [rocket sweep](rocket-parity/comparison.json)
and [Python results](nasa-python-api-results.json) contain the numerical evidence.
Relative differences use nonzero native values; zero fields use absolute bounds.
NaN locations and status codes must match. Some older sweep assertions use
`max(abs(native), 1 SI unit)` scaling; the maximum above was recomputed from
the recorded unscaled values. Trace species and finite-difference noise use
the explicit bounds in their tests.

This establishes parity for the compared inputs. It does not establish
bit-for-bit agreement for every chemistry or compiler.

## Control-flow review

The exhausted-iteration fallback follows `equilibrium.f90` lines 2604–2654,
including its one-use gas reset, enthalpy guard and seed-restoration order.
No reachable input was found for that block: unconverged Newton iterations
cycle before it, while every condensed change clears convergence and resets
`iter` to -1. It must not run on every exhausted Newton epoch. The added
positive-temperature/gas-mole guards after nine singular matrices also have
no isolated native reproducer. The existing water-transition snapshots
already cover the return after too many condensed-set convergence passes;
that return was already leaving the exposed report untouched.

Two proposed corrections would change native behavior and were skipped:

- The rocket pressure-schedule flux uses sonic velocity in NASA's
  `rocket.f90` lines 661–662 and 723–724. Its throat and area solves use flow
  velocity, as Rust already does. Replacing the pressure-schedule basis
  breaks the native FAC early-freeze snapshot. The pressure-schedule area
  comparison now also uses a tighter tolerance.
- NASA's Python `EqSolution` creates a default Fortran solution before
  applying `T_init` and `nj_init` through setters (`CEA.pyx` lines 1546–1563).
  The C setters do not save a seed (`bindc.F90` lines 2859–2871 and
  2960–2989), so the first unconverged solve restores the default seed
  (`equilibrium.f90` lines 2474–2480). Both the successful post-solve state
  and a failed HP state are now checked against native. Passing those
  Python guesses into the first Rust solve would remove this native quirk.

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

### Named gap: singular-ion-shock-roundoff

N2/CH4 at mole fractions 0.985/0.015, 291 K, 0.224 mmHg and 1200 m/s,
with ions, transport and trace 1e-10, has no stable numerical oracle at the
usual 1e-4 relative tolerance. Both solvers return `converged=False` and
`last_error=9`. Their retained station fields remain numerically different.

| Field | Native wheel | Rust before this fix | Rust after weight-scale fix |
|---|---:|---:|---:|
| P21 | 187.1000520352476 | 93.93591779797539 | 93.13409970167557 |
| T21 | 7.4476985618769715 | 6.260995578175356 | 5.854726378439977 |
| T2, K | 2216.705271501316 | 1821.9497132490285 | 1663.7384500747291 |
| Status | 9 | 9 | 9 |

The original P21 difference was 49.7938%; the current difference is 50.2223%.
The changes follow NASA's iteration and initialization arithmetic; they do not bound
the error in this unstable final state.

The Fortran/Rust trace found two defects in failed-shock finalization.
Rust restored the disabled electron equation when rebuilding thermodynamic
partials, and carried `n - sum(nj)` into a diagonal that NASA leaves zero.
The fresh partials matrix now omits that row and leaves that diagonal zero.
The converged Newton-matrix path retains its existing residual arithmetic;
this change concerns the fresh partials requested by shock finalization.
The shock equations now also use the fused products generated by the native
compiler, including the determinant, its two numerators and `1-rho12**2`.
The first TP partials have an ordinary Rust numerical regression:
`dlnV/dlnP=-530396717550.6311`, `dlnV/dlnT=0.9953529048453027`,
`Cp_eq=410380578944840.75 J/(kg K)` and `n=19742862968.172836 kmol/kg`.

`correct_singular` follows the same recovery order, including disabling ions
on an electron-row singularity and reducing element equations only after
ions are inactive. The HP initialization and first nine shock TP iterates
then follow the same path, with small floating-point differences. Neither
an extra recovery solve nor a different seed is needed.

At shock iterate 10, the native determinant `a*d-b*c` is
`3.306607940518319e31`, just `5.3274819966200714e-17` of `abs(a*d)`.
Increasing `a` from `1.125106317305431e19` to its next representable value,
`1.1251063173054312e19`, changes the fused determinant to
`-7.991275834137291e31`. That reverses the pressure/temperature update.
This is cancellation in finite, initialized matrix entries; the observed
difference does not require an uninitialized-memory explanation.

The same NASA source and inputs give these results with GNU Fortran 16.2.0
on macOS arm64. Equilibrium and transport objects were identical in all runs.

| Native shock module | P21 | T21 | T2, K | Status |
|---|---:|---:|---:|---:|
| `-O3`, agrees with the 3.3.4 wheel | 187.1000520352476 | 7.4476985618769715 | 2216.705271501316 | 9 |
| `-O3 -ffp-contract=off` | 43.063163280464536 | 3.9349670429981782 | 1145.0754095124698 | 9 |
| `-O3`, one-ULP change to `G(1,1)` at iterate 10 | 87.79325751251535 | 5.68948239101377 | 1695.25967991425 | 9 |

The default native run exhausts all 60 shock iterates. It retains the last
TP station with the updated P21/T21 ratios; T2 need not equal T1*T21. The
contraction-disabled run instead produces a nonfinite sixth iterate and
restores the fifth TP state and its ratios. Rust's last-valid restoration
already follows that ordering. The comparison checks both T2 and T21,
without imposing an identity that NASA does not preserve on these exits.

These experiments rule out a 1e-4 bound across NASA's own arithmetic choices
for this failed trajectory. Reproducing a particular compiler's last bits
through the entire runaway iteration is still open. The port does not tune
steps, replace species amounts or special-case this mixture to force a final
number. The full numerical comparison remains an explicit expected failure;
an unexpected success fails the suite so the gap must be reviewed.
A wider N2/CH4 sweep from 600 to 2000 m/s also retains status 9 with numerical
differences; the 3000 m/s case converges and matches within 3.24e-10.

[Recorded evidence](singular-ion-shock-roundoff.json) includes the matrix
entries, native compiler experiments, before/after station fields and speed
sweep. The [native-only probe](../scripts/probe_singular_ion_shock.f90)
asserts the last-valid status and positive retained temperature/moles.
Given a release build `B` of the pinned NASA source and a scratch directory
`S`, compile it without changing either repository:

```bash
mkdir -p "$S/probe-modules"
native_includes=(-I"$B/source/modules" -I"$B/extern/fbasics/source/modules")
native_probe() {
  gfortran -O3 "${native_includes[@]}" scripts/probe_singular_ion_shock.f90 "$@" \
    "$B/source/libcea_core.a" "$B/extern/fbasics/source/libfbasics_core.a" \
    "$B/source/bind/c/libcea_abort_support.a" -o "$S/ion-probe"
  "$S/ion-probe" "$B"
}
native_probe
# N is the pinned NASA source directory. Override only the shock object.
gfortran -O3 -ffp-contract=off "${native_includes[@]}" -J"$S/probe-modules" \
  -c "$N/source/shock.f90" -o "$S/shock-off.o"
native_probe "$S/shock-off.o"
```

For the one-ULP run, copy `shock.f90` into scratch, insert
`if (i == 10) G(1,1) = nearest(G(1,1), 1.d0)` immediately before the
incident-equilibrium `X(3)` assignment at source line 518, and compile that
copy with `-O3` and the same module include paths. Pass its object to
`native_probe` in the same position as `shock-off.o`.

EqSolution initial guesses are exposed, including NASA's Python setter
quirk: the first solve restores the default saved seed. NASA's broken
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
