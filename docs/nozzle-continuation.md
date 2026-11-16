# Nozzle continuation checks

The nozzle keeps the previous verified temperature, gas log amounts, stored
amounts, total gas moles, condensed amounts, active phase order and coexistence
pair. NASA resets Newton history and the main trace cutoff on each solve;
the relaxed trace and ion thresholds survive singular recovery. The private
core entry point follows these rules and checks carried phases before assembly.
Public Rust and Python solve signatures are unchanged.

Reference: NASA CEA commit `2f79a647e85737742fb08a127da00d1ea1d069da`,
`source/rocket.f90:RocketSolver_set_init_state` and
`source/equilibrium.f90:EqSolver_solve`. The source was read in a separate clone;
no scripts from that clone were executed. Numerical checks use PyPI `cea==3.3.4`
and its packaged data, independently of Rust. Test comments record the native
API calls and inputs for each golden.

- Air/Jet-A(L), O/F=25, reactants at 298.15 K, 0.7 MPa, area=200:
  chamber and exit match the native rocket within 1e-5 relative. Native stops
  its throat at Mach=0.9999816; Rust preserves its existing sonic closure.
  The raw native throat pressure differs by 2.08e-5 relative and Isp by
  1.62e-5. Refining the native throat with EqSolver SP and Eq. 6.17 until
  `abs(Mach²-1)<1e-10` gives agreement within 1e-5 at every station.
- The same mixture at 20 MPa, area=80 now converges with ice and liquid water.
  Phase adjacency uses the ordered phase records, rather than NASA's compiled
  numbers, which count polynomial intervals too. The independent rerun and
  conservation/phase-equilibrium checks are recorded below.
- RP-1311 example 13 uses N2H4(L)/Be(a)/H2O2(L), weights
  [0.536,0.134,0.33], 298.15 K, 3000 psia, insert BeO(L), trace=1e-10,
  and pressure ratios [3,10,30,300]. All six stations match native at 1e-5
  for T, P, Isp, area, gamma, equilibrium Cp, density and vacuum Isp.
  The explicit trace setting also lets the original checkout solve this deck;
  the older failure report omitted it. Alumina coexistence regressions remain.

The native throat stopping rule and the invalid native 20 MPa exit prevent a
direct 1e-5 comparison against both unmodified native rocket outputs. The
20 MPa regression checks conservation and phase equilibrium instead. No
convergence threshold or existing test tolerance was changed.

A release-build timing on macOS arm64 averaged 100 H2/O2 equilibrium rockets
at weights [1,6], 298.15 K, 3 MPa and area=80. Mean runtime fell from
494.9 µs before continuation to 380.3 µs after it, about 23% faster. This is a
local timing, not a general performance claim.

## 20 MPa water-phase audit

Independently rerun on 2026-11-16 using PyPI `cea==3.3.4`, numpy 2.5.5,
Python 3.15.0 and Rust 1.99.1 on macOS arm64. Exact inputs: Air/Jet-A(L),
mass weights `[25,1]`, inlet temperatures `[298.15,298.15]` K, IAC equilibrium,
chamber pressure 20000000 Pa, supersonic area ratio 80, complete products,
condensed products enabled, ions and transport disabled, default trace and
no inserted/omitted species. Reactant enthalpy is -73913.566459964 J/kg.

At `bc6c53ecbcb2bc04fd77a96d1a7a008fc756f248`, Rust returns
`not_converged: Gibbs/Newton solver failed after 49 iterations`. With
`CEA_DEBUG=1`, it repeatedly replaces liquid at 272.995784 K with ice, then
ice at 273.345644 K with liquid. Ice's compiled phase number is 1; liquid's
is 3 because its record contains two fit intervals. They are adjacent physical
phases. The correction uses their order among all sibling phase records,
including siblings absent from a restricted product pool. The existing
coexistence Newton equations then solve both amounts and temperature.

Native was called independently:

```python
cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
reactants = cea.Mixture(["Air", "Jet-A(L)"])
products = cea.Mixture(["Air", "Jet-A(L)"], products_from_reactants=True)
solver = cea.RocketSolver(products, reactants=reactants, transport=False)
solution = cea.RocketSolution(solver)
weights = np.array([25., 1.])
hc = reactants.calc_property(cea.ENTHALPY, weights, np.array([298.15]*2)) / cea.R
solver.solve(solution, weights, 200., supar=[80.], hc=hc)
```

It prints three `Convergence failed to establish set of condensed species`
warnings, then reports `converged=True`, `last_error=0`. Its reported exit
fields are internally inconsistent:

| Exit property | Native reported | Corrected Rust |
|---|---:|---:|
| Temperature, K | 273.412284434358 | 273.121944434780 |
| Pressure, Pa | 9996.450335475 | 9971.217843957 |
| Density, kg/m³ | 0.152916238648 | 0.129322587823 |
| Enthalpy, J/kg | -1766556.650409 | -1781207.855573 |
| Entropy, J/(kg K) | 7418.561422304 | 7418.561423313 |
| Velocity / Isp, m/s | 1847.754263294 | 1847.860586960 |
| Vacuum Isp, m/s | 1889.585684770 | 1889.586381525 |
| Mach | 6.255539535 | 6.868313044 |
| Sound speed, m/s | 295.378872898 | 269.041404357 |
| Reported area ratio | 80.000062888 | 80.000000033 |
| H2O gas mass fraction | 0.0413212554221 | 0.0374460287379 |
| H2O(L) mass fraction | 0.00630439953709 | 0.00763772631823 |
| H2O(cr) mass fraction | 0 | 0.00254189992170 |

The checks below recompute properties from each solver's reported mass
fractions and the nine-coefficient polynomials in the unchanged
`vendor/thermo-species/data/thermo.inp`. They do not trust the solver's reported
H, S or density. Let `n_j=y_j/M_j` in kmol/kg and let `n_g` sum gas only:

```text
elements_k = sum_j a_kj n_j
h = R sum_j n_j (H_j/R)
s = R [sum_gas n_j (S_j/R - ln(n_j P/(n_g 1e5))) + sum_condensed n_j S_j/R]
rho = P/(n_g R T)
p_H2O = P n_H2O/n_g
p_sat,c = 1e5 exp((G_c - G_H2O,gas)/(RT))
```

Here `R=8314.51 J/(kmol K)`. Air's element counts use NASA's compiler values
N=1.56168, O=0.419590, Ar=0.009365, C=0.000319. The script first checks its
independent H/S sums against the well-converged hot native chamber. The chamber
entropy is 7418.561422364 J/(kg K) in native and 7418.561422207 in Rust;
the chamber stagnation enthalpies are -73913.481210 and -73913.481152 J/kg.

| Independent exit check | Native | Corrected Rust |
|---|---:|---:|
| N element relative residual | -1.17736e-9 | -1.27753e-9 |
| O element relative residual | -1.15163e-9 | -8.02831e-10 |
| Ar element relative residual | -1.18135e-9 | -1.32530e-9 |
| C element relative residual | -1.18138e-9 | -1.32532e-9 |
| H element relative residual | -1.18134e-9 | -8.85742e-10 |
| Recomputed entropy, J/(kg K) | 7457.461482519 | 7418.561413574 |
| Entropy minus chamber, J/(kg K) | +38.900060154 | -8.63318e-6 |
| Recomputed h, J/kg | -1770364.225514 | -1781207.853207 |
| h + u²/2 minus chamber h, J/kg | +10647.164457 | +0.002365865 |
| Recomputed density, kg/m³ | 0.128696817558 | 0.129322587985 |
| Throat mass flux, kg/(m² s) | 19117.609043227 | 19117.609051307 |
| Exit mass flux, kg/(m² s) | 237.800093316 | 238.970113342 |
| Area ratio from these mass fluxes | 80.393614555 | 79.999999933 |
| Relative mass-flux residual at area 80 | -0.00489609228 | +8.38259e-10 |
| Water gas partial pressure, Pa | 671.050157385 | 610.425113311 |
| Saturation pressure over liquid, Pa | 623.420270861 | 610.425113468 |
| Saturation pressure over ice, Pa | 625.176089661 | 610.425113355 |
| G_liquid - G_ice, J/mol | -6.393557801 | +4.20570e-7 |

Even using native's reported properties directly gives an energy residual of
+14454.739562 J/kg and a relative mass-flux residual of +0.182372257. Its
reported area ratio is therefore insufficient evidence of geometric closure.

Liquid is the stable condensed phase at native's reported temperature: it has
the lower Gibbs energy. However, the gas is supersaturated by 7.64009% relative
to that liquid, so the composition is not in phase equilibrium. Rust has
positive ice and liquid amounts, equal Gibbs energies, and vapour pressure
equal to both saturation pressures within 2.6e-10 relative. This satisfies
the water phase-coexistence conditions. The polynomial Gibbs crossing is
273.121944 K, about 0.028 K below the nominal 273.15 K fit boundary, within
NASA's existing 50 K condensed-phase guard band. No new extrapolation rule or
real-water pressure correction was introduced.

A fresh native `EqSolver(products, reactants=reactants)` SP solve with
`state1=solution.entropy[-1]*1000/cea.R`,
`state2=solution.P[-1]` (bar), and the same weights returns
`converged=False`, `last_error=8` (`CEA_NOT_CONVERGED`). Using the independently
recomputed entropy 7457.461482519 instead converges at 274.313008222 K, which
also fails to reproduce the rocket's reported exit. Fresh Rust SP solves
at its chamber entropy converge at both pressures: at native's 9996.450335475 Pa
they give 273.121944437 K; at Rust's 9971.217843957 Pa they give the same
temperature and reproduce all three water amounts within 1.6e-10 mass fraction.

The evidence supports a sound equilibrium state in this model, and rejects
native 3.3.4's reported success as a usable golden for this case. The regression
`air_jet_a_high_pressure_exit_obeys_conservation_and_water_phase_equilibrium`
checks the independent balances, both chemical potentials, and a fresh SP
solve. It fails with the original core's 49-iteration error and passes with
the correction.

The before/after sweep covers O/F=25 at 298.15 K and O/F=15 at both 298.15 and
300 K, pressures [0.7,3,7,10,20] MPa, and areas [15,20,50,80,100,200]: 90 cases.
All 89 previously successful cases have bit-identical state and flow fields
and compositions; only the target error becomes a success. SHA-256 fingerprints
of all stations are emitted by the audit script for a repeatable comparison.

Reproduce from the repository root, without changing `mise.toml`:

```sh
S=$(mktemp -d)
mise exec -- env UV_CACHE_DIR="$S/uv-cache" TMPDIR="$S" uv venv "$S/venv"
mise exec -- env UV_CACHE_DIR="$S/uv-cache" TMPDIR="$S" \
  uv pip install --python "$S/venv/bin/python" cea==3.3.4 numpy maturin
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  "$S/venv/bin/maturin" build --release --out "$S/wheels"
mise exec -- env UV_CACHE_DIR="$S/uv-cache" TMPDIR="$S" \
  uv pip install --python "$S/venv/bin/python" "$S"/wheels/*.whl
"$S/venv/bin/python" scripts/check_air_jet_a_nozzle.py --output "$S/after.json"
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo test --release
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo fmt --check
```

To reproduce the baseline, install a wheel built from `bc6c53e` into the same
throwaway venv and run the same audit script to `before.json`. It records the
error rather than hiding it. Compare corresponding `neighbors` objects:
89 pairs are identical and the target changes from `error` to
`stations_sha256`. This audit used the same interpreter and dependencies on
both sides.

Validation: `cargo test --release` passes 58 tests with one pre-existing ignored
diagnostic and no failures; `cargo clippy --release --all-targets -- -D warnings`,
`cargo fmt --check`, and Ruff on the audit script pass without warnings.
