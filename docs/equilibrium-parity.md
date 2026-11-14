# Equilibrium parity evidence

This branch matches NASA CEA 3.3.4 for the cases below. These results supersede
the older temperature-floor, N2O-reference and omitted-condensed-record notes
in the README and implementation matrix's historical prose.

The implementation uses NASA's 160–22000 K equilibrium bounds. Gas fits are
selected by their lower bounds and extrapolated as NASA does. Two-interval
product gases receive NASA's generated high-temperature extension. Nonzero
extrapolated products are listed in `State.extrapolated_species`; the existing
`trace_extrapolated_species` remains the subset at mole fraction <= 1e-11.
The fit-boundary comparison tolerance is 1e-3 K.

Fitted reactant properties extrapolate. Reference-only `Reactant.properties`
(`thermo` in Python) accepts NASA's reported reference range, T_ref ± 10 K,
and returns no Cp or entropy. Stream `reactant_enthalpy` matches NASA's
`mixture_calc_enthalpy`: it uses the stored enthalpy for a reference-only
reactant at any finite positive inlet temperature. No heat capacity is invented.
Air uses the precise element counts restored by NASA's database compiler.

Condensed phases use NASA's compiled interval/phase numbering. This fixes
solid/liquid alumina oscillation. At coexistence, Cp_eq is zero by NASA's
reporting convention; gamma is calculated from the pressure derivative after
removing the liquid phase from the derivative matrix. All 11 previously
omitted species are retained, including their reversed/zero-width intervals.
Duplicate-name intervals remain in source order; their bounds use min/max
across the intervals, as in NASA.

Rust adds `Solver::with_omit` and `Solver::with_trace`. Python accepts `omit`
and `trace` constructor keywords. Omit uses NASA's exact-name filtering; names
absent from the database have no effect. An omitted inserted phase and invalid
trace values fail explicitly. Smooth truncation is omitted:
NASA defaults to hard truncation, and all reported comparisons pass with it.

The source checkout was cloned outside this repository and inspected at
`2f79a647e85737742fb08a127da00d1ea1d069da`. No reference scripts were executed.
The numerical oracle was the separately installed native PyPI `cea==3.3.4`.

Validation on macOS arm64, Rust 1.99.1 and Python 3.15:

- 111 recorded states compare T, P, density, H, S, Cp_eq, gamma_s and M at
  1e-4 relative tolerance (absolute tolerance for a zero reference Cp), and
  mole fractions at 3e-5 absolute tolerance. Worst differences were
  1.571e-8 relative and 6.954e-9 absolute respectively.
- These states include AP/aluminium at mass O/F 4 and 6 and AP/aluminium/
  butadiene at 68/18/14, each at 0.7, 3, 7 and 20 MPa; air/Jet-A at O/F 15
  at 7 and 20 MPa with both 298.15 and 300 K inlets; and LOX/LH2 at O/F 3 at 3 MPa. Each has an HP chamber
  and independent SP solves at native nozzle exit pressures for areas
  15, 50 and 200. They also include every equilibrium point in RP-1311
  examples 1, 2, 3, 4, 5 and 14, added condensation points at 200, 180 and
  160 K, and H2/O2 and CH4/O2/N2 at 6000.1, 10000, 20000 and 22000 K.
- A separate 135-nozzle sweep compares all three stations directly with
  native RocketSolver. AP/aluminium O/F 4 and 6 and the 68/18/14 blend use
  0.7, 1, 3, 5, 7, 10 and 20 MPa; air/Jet-A uses 7, 10 and 20 MPa;
  LOX/LH2 O/F 3 uses 0.7, 3 and 20 MPa. Every mixture/pressure uses areas
  15, 20, 50, 100 and 200. All converge in both implementations. Worst
  property difference is 5.603e-5 relative (Jet-A exit pressure); worst
  mole-fraction difference is 5.462e-6 absolute (alumina coexistence).
- The original `scripts/compare_native.py` passes all 146 rows. Worst scaled
  difference, using that script's max(abs(native), 1 SI unit) scale, is
  5.841e-5.
- `cargo +1.99.1 test --locked` and
  `cargo +1.99.1 clippy --features python -- -D warnings` pass. The existing
  ignored upstream diagnostic sweep remains ignored.
- The wheel tests' three other methods pass. `test_install_and_bad_inputs`
  fails at its old assertion that N2O at 298.15 K raises DomainError.
  That assertion must be changed by the owner of `tests/test_wheel.py`.

Regenerate the committed fixtures in a throwaway uv project with
`cea==3.3.4 numpy` installed:

```sh
python scripts/record_equilibrium_parity.py --nasa-source /tmp/nasa-cea
cargo +1.99.1 test --locked native_equilibrium_domain_regressions -- --nocapture
```

The generator only reads sample literals with `ast.literal_eval`. Historical
sample omission lists contain misspelled or reference-only names which NASA
ignores; the fixtures preserve those lists as written. Other
reactant/property goldens in `tests/solver.rs` document their native API calls.

Integration hooks:

- In `rocket.rs`, change `.fold(300.0_f64, f64::max)` to
  `.fold(160.0_f64, f64::max)` in the frozen entropy root's `min_t` calculation.
  This branch leaves that line to the rocket owner. The only edits in that
  file refresh the new extrapolation metadata in cloned frozen states.
- Frozen expansion still applies the existing per-species range rejection;
  matching NASA's gas extrapolation and condensed frozen-phase limits belongs
  with that rocket change. The 135-nozzle sweep above uses equilibrium flow.
- Update the wheel test's obsolete N2O rejection and the README's historical
  domain statements during integration. README was deliberately not edited.
