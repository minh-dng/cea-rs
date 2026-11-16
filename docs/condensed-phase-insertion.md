# Condensed-phase insertion checks

Compared with `ff4c575` on 2026-11-16–17, using Rust 1.99.1, Python 3.15.0,
numpy 2.5.5 and independent native `cea==3.3.4` on macOS arm64. All cases use
100000 Pa, complete products, condensed products enabled, ions disabled,
default trace, and no insert/omit list unless stated. No source data or
convergence tolerance changed.

## Causes and corrections

1. A cold HP target can lie below every gas-only state, so the Newton loop
   never converges and its gated condensed insertion test never runs. SP can
   instead find an invalid gas-only root below 160 K. After an unsuccessful
   gas-only solve, obtain a converged TP state at NASA's 160 K lower bound
   (TV for volume constraints), then retry the original constraints with its
   condensed active set. Retry only once, and only when that seed converges
   with a condensate. Existing converged solves keep their original path.
2. Br2(cr) declares 300 to 265.9 K. Requiring both bounds excludes every
   temperature. Keep the source interval and public range metadata, use its
   upper phase limit, and disregard the unusable lower bound for insertion.
   Swapping endpoints would still exclude 262 K. This uses the unchanged
   polynomial by extrapolation; the State still reports that extrapolation.
3. Sibling condensates have identical element columns. Without a temperature
   unknown their joint Newton matrix is singular. For TP/TV, choose the
   lower-Gibbs sibling, transferring the amount instead of creating a
   coexistence pair. A sibling is compared only inside its own fit range:
   extrapolated fits can look falsely stable. Without that limit NaCN(III)
   replaced the stable NaCN(II) at 240 to 288 K, and Cr2O3(I) replaced
   Cr2O3(I') at 304 K; a scan of the database found ten such sibling pairs.
   Equal Gibbs energies keep one phase; a fixed-T phase split is
   undetermined. HP/SP phase coexistence and its temperature equation stay
   unchanged.

   Fixed-T continuation from the other phase and invalid inserted phases are
   covered by the recovery audit below. No public entry point starts a fixed-T
   solve from a state; rocket continuation uses SP.

## NASA source and native behavior

Reference: NASA CEA commit `2f79a647e85737742fb08a127da00d1ea1d069da`.
Only the downloaded source text was inspected; no reference scripts were run.

- [EqSolver_solve](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/equilibrium.f90#L2541)
  cycles before `test_condensed` when the base Newton solve has not converged.
  Native unseeded Hg HP returns `converged=False`, error 8. The native run
  with `insert=["Hg(L)"]` also returns error 8; Rust's already-working seeded
  HP/SP results retain every field and composition bit for bit.
- [build_thermo_coeffs](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/thermo.f90#L604)
  copies the reversed interval without normalization.
  [get_valid_temperature_range](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/thermo.f90#L117)
  takes the minimum of the first column and maximum of the second;
  [insertion](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/equilibrium.f90#L1913)
  uses those columns directly. Native Br2 TP reports success with no solid,
  h=175585.982875035 J/kg and p_Br2=83331.942659745 Pa. That is
  18.2449054 times solid saturation (16.9947061 times liquid saturation).
  The Rust change deliberately corrects this impossible insertion condition.
- [check_condensed_phases](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/equilibrium.f90#L1648)
  runs for fixed-T solves too, adds sibling phases and assigns the nominal
  transition temperature without a fixed-T or Gibbs-ordering guard. Native
  TP fails with error 8 for water at 273.15 K, Hg at 234.29 K and NaCN at
  288.5 K. At 835 K native NaCN succeeds with liquid, liquid mass fraction
  0.4998693263144834 and p_NaCN/p_sat,L=0.999999999999988.

Native inputs use fresh `Mixture(["Ar", name])`, an independent
`Mixture(["Ar", name], products_from_reactants=True)`, `EqSolver`, and
`EqSolution` for each row. Calls use `solve(sol, cea.HP, h/cea.R, 1., weights)`
or `solve(sol, cea.TP, T, 1., weights)`; results are checked with `converged`
and `last_error`, not with returned property values alone.

## Independent first-principles checks

Read coefficients directly from unchanged `thermo.inp` with the independent
parser in `scripts/check_air_jet_a_nozzle.py`. Evaluate the nine-coefficient
polynomials, rather than using the reported h or s. With `n_j=y_j/M_j`
(kmol/kg), gas-only `n_g`, and R=8314.51 J/(kmol K):

```text
b_k = sum_j a_kj n_j
h = sum_j n_j H_j
s = sum_gas n_j [S_j - R ln(n_j P/(n_g 1e5))] + sum_condensed n_j S_j
G_j = H_j - T S_j
p_sat,c = 1e5 exp((G_c - G_vapor)/(R T))
mu_vapor/(R T) = G_vapor/(R T) + ln(p_vapor/1e5)
```

Ar/Hg and Ar/Br2 weights are [1,20]; the four transition cases use [1,1].
The HP target is -15551.803948663037 J/kg. The SP target is
484.1589263235252 J/(kg K), the mean of the 232 K and 237 K TP entropies.

| Case | T, K | Positive condensed mass fractions | Recomputed assigned quantity minus target |
|---|---:|---|---:|
| Hg HP | 234.28999933688 | Hg(L): 0.479137243881, Hg(cr): 0.473243707791 | +5.82077e-11 J/kg |
| Hg SP | 234.28999933688 | Hg(L): 0.478949938993, Hg(cr): 0.47343101268 | +3.18323e-12 J/(kg K) |
| Br2 TP | 262 | Br2(cr): 0.943263839016 | +0 K |
| Water TP | 273.15 | H2O(L): 0.498612290857 | +0 K |
| Hg TP | 234.29 | Hg(L): 0.499999992559 | +0 K |
| NaCN TP (288.5 K) | 288.5 | NaCN(III): 0.499999999526 | +0 K |
| NaCN TP (835 K) | 835 | C(gr): 3.14031995969e-05, NaCN(L): 0.499869326341 | +0 K |

Element residuals below are signed `(product/feed)-1`, recomputed from the
normalized reported mass fractions. Minor differences from `element_residual`
come from that normalization; none changes a feed or convergence check.

| Case | Element relative residuals |
|---|---|
| Hg HP | AR: -7.77156e-16, HG: +0 |
| Hg SP | AR: -1.66533e-15, HG: +0 |
| Br2 TP | AR: -2.77556e-15, BR: +2.22045e-16 |
| Water TP | AR: -3.12945e-09, H: -3.12946e-09, O: +3.91806e-09 |
| Hg TP | AR: +4.44089e-16, HG: -3.33067e-16 |
| NaCN TP (288.5 K) | AR: -9.47902e-10, NA: -9.47902e-10, C: +3.48171e-09, N: +1.88683e-09 |
| NaCN TP (835 K) | AR: -3.59859e-11, NA: +8.78098e-12, C: -3.59859e-11, N: +1.42353e-10 |

For coexistence both phases must have equal chemical potential. For a single
phase, the absent sibling must have higher G, and vapour must be saturated
with the present phase and undersaturated with the absent sibling.

| Case | G_second - G_first, J/mol | p_vapour, Pa | p/p_sat (first; second) |
|---|---:|---:|---|
| Hg HP | +2.98023223877e-11 | 0.00029637996232826 | Hg(cr): 1.00000000000015; Hg(L): 1.00000000000013 |
| Hg SP | +2.98023223877e-11 | 0.000296379962328307 | Hg(cr): 1.00000000000031; Hg(L): 1.0000000000003 |
| Br2 TP | +154.631654277 | 4567.40886333138 | Br2(cr): 1.00000000000086; Br2(L): 0.931476798338866 |
| Water TP | -0.617297710478 | 611.670267367096 | H2O(cr): 0.999728224743234; H2O(L): 0.999999992404983 |
| Hg TP | -6.49644434452e-06 | 0.00029637998896418 | Hg(cr): 0.999999996665474; Hg(L): 1.0000000000004 |
| NaCN TP (288.5 K) | -31.6390849032 | 1.05043151607087e-22 | NaCN(II): 0.98689670292013; NaCN(III): 0.999999995325958 |
| NaCN TP (835 K) | -1.35313570499e-05 | 0.413864818697597 | NaCN(III): 0.99999999794196; NaCN(L): 0.999999999890987 |

The sibling order is Hg(cr)/Hg(L), Br2(cr)/Br2(L), H2O(cr)/H2O(L),
NaCN(II)/NaCN(III), or NaCN(III)/NaCN(L). Hg coexistence has
G_L-G_cr=2.98023e-11 J/mol and maximum |(mu_c-mu_v)/(RT)|=3.091e-13.
NaCN's reported vapour is tiny at 288.5 K but positive:
p_NaCN=1.05043151607e-22 Pa. Its saturation ratio was checked independently
as well.

`tests/condensed.rs` contains one regression per cause, with shared composition
sums. The transition regression also checks TV, the exact water Gibbs crossing
at 273.121944434780 K, and the exact Hg crossing at 234.28999933688027 K.
The test enforces Gibbs ordering tightly enough to distinguish the 6.50e-6
J/mol Hg gap and the 1.35e-5 J/mol NaCN gap. All three final tests were run
against the original `ff4c575` core: 0 passed, 3 failed. With the corrected
core: 3 passed, 0 failed.

## Unchanged results and checks

- All 90 Air/Jet-A neighbors in `scripts/check_air_jet_a_nozzle.py` succeed
  before and after, with identical SHA-256 station fingerprints. The 20 MPa,
  area-80 coexistence exit and standalone SP rows also compare exactly.
- All 120 H2/O2 nozzle runs compare bit for bit: H2/O2 at weights [1,6],
  298.15 K, and H2(L)/O2(L) at weights [1,5.55157], 20.27/90.17 K;
  pressures [0.7,3,5.33172,20] MPa, areas [1,1.001,10,80,200], and
  equilibrium/frozen-at-chamber/frozen-at-throat flow. Compare every station's
  T, P, rho, h, u, s, both Cp values, gamma, molar mass, both compositions,
  velocity, sound speed, Mach, area, c-star, thrust coefficient and vacuum Isp.
- All 111 equilibrium fixture inputs compare exactly, including the above
  state properties, iterations, element/constraint residuals and extrapolation
  lists. Seeded Hg HP and SP also keep all fields and compositions unchanged.
- No changed number was found in those existing results. No fixture was updated.
- Baseline Rust suite: 58 passed, 0 failed, 1 ignored. Final release suite:
  61 passed, 0 failed, 1 ignored. The existing upstream diagnostic sweep stays ignored.
- Installed-wheel suite: 4 passed before and after. The release Python wheel
  was rebuilt and tested outside the source package.
- `cargo clippy --release --all-targets -- -D warnings`: passed, no warnings.
- `cargo fmt --check`: passed.

Build commands used a local temporary directory; `mise.toml` was unchanged:

```sh
S=$(mktemp -d)
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo test --release
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo fmt --check
```

The verification environment also overrode `UV_CACHE_DIR` for native CEA,
numpy and maturin. Native CEA never supplied a seed or fallback result to Rust.
This audit covers these mixtures and sweeps; it does not establish validity of
extrapolated Br2 fits or broad metal-system parity.

## Fixed-temperature recovery audit against f873de1

Checked on 2026-11-17 with Rust 1.99.1, Python 3.12.15, numpy 2.5.5 and an
independent `cea==3.3.4` wheel. All inputs below use Ar plus the named reactant,
weights [1,1], 100000 Pa, full products, condensed enabled, ions disabled and
default trace. No database coefficient, convergence tolerance or public API
changed.

The original numerical path runs first. Only a failed TP/TV solve retries
once with these corrections; the retry must pass the same Newton and public
residual checks. Its reported iteration count includes both attempts.

1. **Unseeded NaCN at 300 K:** after forming NaCN(III), a tiny negative graphite
   amount is removed. The remaining matrix has a nonzero trace-element pivot,
   but its pivot ratio exceeds the Gaussian solver's `1e25` sentinel. Treating
   that pivot as singular lets G.2 remove the valid NaCN phase. The retry
   accepts finite nonzero pivots without that arbitrary ratio limit, in both
   Newton and the final property partials. Zero/nonfinite pivots still fail.
2. **Inserted NaCN(II) at 300 and 835 K:** the seed reaches Newton before its
   fit is checked. At 300 K its extrapolated Cp is negative, so F.3.1 removes
   it despite its positive amount. At 835 K the far-extrapolated fit prevents
   convergence. F.4 also lets an invalid incumbent's lower extrapolated Gibbs
   energy veto its valid sibling. On the retry, validate inserted phases
   before matrix assembly; an invalid incumbent cannot veto an allowed phase.
   Compare Gibbs energies when both phases are allowed.
3. **Ice at 270 K continued to 274 K:** F.4 correctly transfers ice to liquid,
   but the carried gas is truncated to Ar and H2O. Its H and O rows are
   dependent. G.2 disables truncation and immediately removes/blocklists the
   only condensate, before trying the restored gas. The retry rebuilds once
   without truncation before any phase removal. The liquid then converges.

NASA's [insertion rule](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/equilibrium.f90#L1912)
requires T at or below the upper fit bound and at or above the lower bound,
except that a 200 K lower bound permits lower temperatures. The port's existing
reversed-interval exception is unchanged. NaCN(II) declares 197.7–288.5 K,
NaCN(III) 288.5–835 K, and NaCN(L) 835–6000 K. Thus **only NaCN(III) is allowed
at 300 K**. Independently evaluated G/(RT) is -50.727820314964 for II and
-50.612971708236 for III: II's apparent 286.472966736 J/mol advantage is an
out-of-range extrapolation, not a competing allowed phase. NASA's 50 K phase
switch neighborhood is not an extension of the insertion fit range.

At 835 K, III and liquid are both allowed. Liquid is lower by
0.0000135313570499 J/mol. At 274 K only liquid water is allowed among the water
condensates; its range is 273.15–600 K. Extrapolated ice also has higher G by
19.3714066682 J/mol, but does not enter the allowed-phase comparison.

Independent checks use the raw `thermo.inp` parser and polynomial sums above,
not the solver's reported residuals or its thermo evaluation functions.
Every assigned T and P is retained exactly. Every positive condensate is
inside its own fit: NaCN(III), NaCN(L), H2O(L), and graphite (300–6000 K).

| Case | Iterations, original + retry | Positive condensed mass fractions | p_vapour, Pa | p/p_sat, present NaCN/water phase |
|---|---:|---|---:|---:|
| NaCN(II) or NaCN reactant, unseeded, 300 K | 114 + 55 | NaCN(III): 0.499999999999995 | 2.08162731570701e-21 | 0.999999999999957 |
| NaCN reactant, insert NaCN(II), 300 K | 96 + 32 | NaCN(III): 0.499999994787154; C(gr): 1.38433917840e-25 | 2.08162725381046e-21 | 0.999999970265262 |
| NaCN reactant, insert NaCN(II), 835 K | 75 + 28 | NaCN(L): 0.499869326615283; C(gr): 3.14030871179e-05 | 0.413864818234227 | 0.999999998771369 |
| H2O(L) reactant, 270 K seed → 274 K | 9 + 4 | H2O(L): 0.498523671941702 | 650.478058624024 | 0.999999999426575 |

The absent allowed NaCN(III) at 835 K has p/p_sat=0.999999996822341 and the
positive Gibbs gap given above. At 300 and 274 K there is no absent allowed
sibling. None of these vapours is supersaturated.

| Case | Signed element residuals, (product/feed)-1 |
|---|---|
| Unseeded NaCN, 300 K, either reactant | AR: -9.32587e-15; NA: +3.08642e-14; C: -9.32587e-15; N: -9.32587e-15 |
| Insert NaCN(II), 300 K | AR: -1.04257e-08; NA: -1.04257e-08; C: -1.04257e-08; N: +6.25299e-08 |
| Insert NaCN(II), 835 K | AR: -4.05665e-10; NA: +9.79736e-11; C: -4.05666e-10; N: +1.60641e-09 |
| Water continuation, 274 K | AR: -1.62315e-10; H: -1.62316e-10; O: +2.03218e-10 |

Native CEA was installed with `uv venv` and `uv pip install cea==3.3.4`, with
`UV_CACHE_DIR` and `TMPDIR` redirected to a temporary directory. Fresh native
reactant/product mixtures, solvers and solutions were used for each input.
All requested native NaCN cases report `converged=True`, `last_error=0`:

| Native case | NaCN condensed mass fraction | p_NaCN, Pa | p/p_sat of that phase |
|---|---:|---:|---:|
| Unseeded, either reactant, 300 K | III: 0.499999999999998 | 2.08162731570703e-21 | 0.999999999999968 |
| Insert II or III, 300 K | III: 0.499999961481994 | 2.08162685834947e-21 | 0.999999780288418 |
| Insert II or III, 835 K | liquid: 0.499869326314486 | 0.413864818742703 | 0.999999999999975 |

Native's largest independent element residual in these rows is +4.62037e-07
for nitrogen in seeded 300 K runs. Native output was never a seed or fallback
for the Rust solver.

### Regression and unchanged-output checks

- One regression per failure: two public tests in `tests/condensed.rs` and one
  crate-internal continuation test in `src/solver.rs`. All three fail at the
  original solve on `f873de1`; the seeded 835 K failure was also reproduced
  independently (75 iterations). The tests check the assigned T, phase,
  balances and vapour chemical potential.
- All 358 previously successful complete snapshots are bit for bit identical:
  111 equilibrium fixture inputs, 120 H2/O2 nozzle runs described above, and
  127 TP phase/seed boundary cases for NaCN, water, Hg and Br2. Compare all
  state fields, both compositions, iteration counts, residuals, extrapolation
  lists and nozzle station fields. No fixture changed.
  The canonical JSON SHA-256 of these 358 snapshots is identical before/after:
  `c16b1570848d856e23fb70fb5d892a40872a1292f964f4fbf5dce88be71dff36`.
- All 90 Air/Jet-A fingerprints in `scripts/check_air_jet_a_nozzle.py` remain
  identical. Its entire JSON report, including the phase exit and standalone
  SP rows, is identical too.
  Report SHA-256:
  `a103f237e2187b43ea88856a227de63e24f71e2549fe331aca1775f5d8fcd0c7`.
- Of 17 failing boundary controls, 12 NaCN cases now converge: 300 K with no
  seed or insert II/III/liquid; 310 K with II; 338.5 K with II/III/liquid; and
  785, 834, 835 and 836 K with II. Every new result passed the independent
  fit, Gibbs, element and no-supersaturation checks. Maximum element residual
  is 4.59784e-07 and maximum allowed-phase saturation ratio is below 1.
- Five unseeded Ar/Cr2O3(I') controls at 300, 304, 305, 320 and 330 K still
  fail. Their error messages now count both attempts: 46→112, 46→112, 50→116,
  49→115 and 49→115 iterations. No successful result elsewhere changed.

The baseline release suite passed 70 tests with 1 ignored. The final release
suite passes 73 tests with 1 ignored. The deck tests pass 3/3 and condensed
tests 5/5; the continuation regression also passes. Both release Clippy runs
(`--all-targets`, default and `--features python`) pass with `-D warnings`, and
`cargo fmt --check` passes. The rebuilt installed-wheel tests pass 8/8 and the
NASA API tests with independent native CEA pass 9/9. Build output and temporary
files use a temporary directory; `mise.toml` is unchanged.

```sh
S=$(mktemp -d)
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo test --release
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo test --release --test nasa_decks --test condensed
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" \
  cargo clippy --release --all-targets --features python -- -D warnings
mise exec -- env CARGO_TARGET_DIR="$S/target" TMPDIR="$S" cargo fmt --check
```
