# NASA input decks and reports

Source: [NASA CEA 2f79a647e85737742fb08a127da00d1ea1d069da](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da),
version 3.3.4, Apache-2.0. `NOTICE` is an unchanged copy of upstream `NOTICE.txt`.
See the repository's root LICENSE, NOTICE and `docs/provenance.md`.

| Local files | Upstream source |
|---|---|
| `example1.inp` … `example14.inp` | `test/main_interface/example*.inp` |
| `example1.out` … `example14.out` | `test/main_interface/reference_output/example*.out` |
| `rp1311_examples.inp` | `samples/rp1311_examples.inp` |
| `sample-example1.inp` | `samples/example1.inp` |
| `native-3.3.4/example*.out` | Unmodified NASA executable, run on the corresponding local `.inp` |
| `native-3.3.4/sample-shock.out` | Combined sample's seventh block, with `output long` added |
| `native-3.3.4/sample-12.out` | Combined sample's twelfth block (additional HNO2 product), with `output long` added |
| `shock-reflected-transport.inp` | Additional regression: example 7 mixture at 1400 m/s, incident/reflected equilibrium/frozen modes and transport |
| `native-3.3.4/shock-reflected-transport.out` | Unmodified pinned NASA executable, run on that regression deck |
| `python-api-3.3.4.json` | Successful calls to PyPI cea==3.3.4 through the Python test helpers; no Rust-generated expected values |

Original decks/reports are copied byte-for-byte, including trailing spaces,
missing final newlines and line endings. `.gitattributes` preserves those bytes.
`SHA256SUMS` covers all `.inp`, `.out` and the NASA notice.
No NASA runtime or network is needed to run the Rust tests.

The Python snapshot was captured on 2026-11-17 with CPython 3.12.14, NumPy
2.5.5 and macOS arm64. Its metadata includes the NASA source-distribution hash.
It contains the 11 successful compared `test_rocket_schedules` inputs through
`rp_rocket(native, ...)`,
all 16 `test_shock_and_detonation` flag combinations through `shock(native, ...)`,
and that test's four `(temperature, pressure)` detonation inputs through
`detonation(native, ...)`. The remaining ten cases are the explicit RP1 TP/RP2
TV inputs and eight RP14 temperatures in
`test_rp_equilibrium_and_condensed_exceptions` in `tests/test_nasa_api.py`.
Each snapshot records the native solution's THERMO fields and application/
transport fields used by that test. Frozen rocket Cv and reflected-frozen
shock Cp/Cv are omitted because NASA retains or leaves those fields
uninitialized; they are never replaced by Rust values. All native calls
converged with `last_error == 0`. Regenerate only by rerunning those same
helpers and inputs against NASA 3.3.4 and serializing its numeric properties.
Both live-native and NASA-free test runs compare these stored values.

The historical reports identify CEA2 (2004). They differ from the modern
executable in a few cells. Example 12's historical report echoes a frozen-only
input; the newer `.inp` requests both equilibrium and frozen branches.
Both sets are retained. The tests record individual known differences and
enforce `1e-4` on every other compared value. See `docs/input-decks.md` for
the passing/failing table, exact differences and comparison policy.

To regenerate modern reports, check out the pinned NASA source outside this
repository. Build with C bindings enabled: this version's Fortran executable
needs its C abort support even when no Python binding is built.

```sh
cmake -S "$NASA_SOURCE" -B "$NASA_BUILD" \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release \
  -DCEA_BUILD_TESTING=OFF -DCEA_ENABLE_BIND_C=ON \
  -DCEA_ENABLE_BIND_CXX=OFF -DCEA_ENABLE_BIND_PYTHON=OFF \
  -DCEA_ENABLE_BIND_MATLAB=OFF -DCEA_ENABLE_BIND_EXCEL=OFF
cmake --build "$NASA_BUILD" -j 4
```

Copy the original 14 decks to a temporary working directory, then run each
there (the `.out` file is written beside its input):

```sh
"$NASA_BUILD/source/cea" -t "$NASA_BUILD/thermo.lib" \
  -r "$NASA_BUILD/trans.lib" example1.inp
```

Repeat for examples 2–14 and copy the reports into `native-3.3.4/`. For the
combined sample's shock report, extract its seventh `end`-terminated block
unchanged, insert `output long` immediately before `end`, name that temporary
input `sample-shock.inp`, and run the same command. Its 1000 m/s equilibrium
column is explicitly marked nonconverged by NASA; the tests exclude only that
column and compare the other 13 shock cases. Repeat this extraction for the
twelfth block as `sample-12.inp`; its additional HNO2 product requires its own
report rather than using the main-interface example 12 report.

Run `shock-reflected-transport.inp` through the same executable to regenerate
its report. It adds the reflected and transport rows absent from the shipped
shock decks. NASA prints equilibrium reflection in a separate section without
repeating the incident state. Its frozen reflected header Cp retains the
incident value; the separate transport Cp uses the reflected temperature.
Both header differences are bounded explicitly in the regression test.

The modern example 6 report is retained for provenance but is not an oracle:
its `UNBURNED GAS` table contains burned properties, followed by a rounded
duplicate labeled `BURNED GAS`. Tests use the complete historical long report
for example 6 instead.

Run the comparisons, including the recorded discrepancies:

```sh
cargo test --release --test nasa_decks -- --nocapture --test-threads=1
```
