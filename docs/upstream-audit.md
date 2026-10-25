# Upstream audit: NASA CEA + existing Rust implementations

Source-audit notes gathered by a free MiMo research task. The task wrote this
file but ended with a free-usage-limit error before returning a completed
result. This is research evidence, not the independent implementation review.
The implementation owner checked the NASA/Rust identities, licenses, format
and solver depth directly. The current implementation matrix is authoritative
for this package; the tables below describe upstream projects before adaptation.

Research only (no code, no commits, no production changes). All clones/caches on T7.
Runs used `CARGO_HOME="/Volumes/T7 Shield/keplr-work/cache/cargo"`.

## 1. Primary sources (local path @ immutable ref)

| Source | Path (T7 `upstream/`) | SHA | Date |
|---|---|---|---|
| NASA CEA https://github.com/nasa/cea | `nasa-cea` | `2f79a647e85737742fb08a127da00d1ea1d069da` (main, shallow) | 2026-09-25 |
| aselimov/cea-rs | `aselimov-cea-rs` | `44e2426aa285a883c6c63937e5de8592475ecbd1` (shallow) | 2026-03-29 |
| bit-mango/cea-rs | `bit-mango-cea-rs` | `914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad` (shallow) | 2026-08-10 |
| bit-mango/thermo-species | `thermo-species` | `d4061dddda7f04e2edf75cdb235d19d1fed4161b` (shallow) | 2026-08-10 |
| mcurrie99/nasa-cea | `mcurrie99-nasa-cea` | `0c85d326c16406f26250fc15992d034d3ed62b7d` (shallow) | 2026-05-19 |

Registry/search evidence (cached `cache/*.json`, fetched 2026-10-25):
`https://crates.io/api/v1/crates/cea-rs` · `/thermo-species` · `/cea` ·
`https://api.github.com/search/repositories?q=cea+nasa+language:rust` (3 hits) and
`...?q=chemical+equilibrium+gibbs+language:rust` (0 hits).
Cited references: NASA TP-2002-211556 `https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf` ·
RP-1311 Pt1 `https://ntrs.nasa.gov/citations/19950013764` · Pt2 `https://ntrs.nasa.gov/citations/19960044559` ·
docs `https://nasa.github.io/cea/`.

## 2. NASA license + AGENTS.md constraints

- Code: **Apache-2.0** (`LICENSE.txt`); `NOTICE.txt` copyright (c) 2025 NASA, Apache-2.0.
- Data attribution in `NOTICE.txt`: `data/thermo.inp`, `data/trans.inp`, `*.lib` → TP-2002-211556;
  `samples/*.inp|*.out` → RP-1311 Part 2. Reuse ⇒ keep `LICENSE` + `NOTICE`, state modifications.
- `AGENTS.md`: do not modify `data/thermo.inp`/`trans.inp`; *"preserve bitwise results and
  scientific behavior"*; avoid algorithmic changes; layer boundaries `source/`, `source/bind/c/`,
  `source/bind/python/`. ⇒ port data byte-for-byte and differential-test against RP-1311 outputs.

## 3. Full mode inventory (`source/input.f90:58`)

Problem types: `tp, hp, sp, tv, uv, sv, det, rkt, shk` (aliases `pt/ph/ps/vt/vu/vs`, `ro`, `sh`).
- Equilibrium: TP/HP/SP/TV/UV/SV (C enum `CEA_TP..CEA_SV`).
- Rocket `rkt`: `eq` vs `fr|fz` frozen, `nfz` frozen iterations, `sub`/`sup` area ratios,
  `pi|pc` pressure, `ac` throat, `fac` finite-area throat, `mach`, `mdot`, schedules `o/f, f/o, f/a,
  phi, %f, h/r, s/r, u/r, rho`.
- Shock `shk`: `inc` incident, `ref` reflected, `mach` (stations 1, 2, 5).
- Detonation `det`: CJ, optional `frozen`.
- Species/solver opts: `ions`, `trace`, `only/omit/insert`, `transport`, `smooth_truncation`
  (default false), `truncation_width` (default 0.25), `tces`.

## 4. API surface

- Fortran core: `equilibrium.f90` 5763 LOC, `main.f90` 4028, `rocket.f90` 1822, `mixture.f90` 1290,
  `shock.f90` 1193, `input.f90` 1031, `thermo.f90` 807, `detonation.f90` 315, `transport.f90` 376,
  `database_compile.f90`, `fits.f90`, `atomic_data.f90`, `units.f90`.
- C API `source/bind/c/cea.h`: **112 distinct `cea_*` functions**; opaque handles; `cea_solver_opts`
  `{trace, ions, transport, reactants, ninsert, insert, smooth_truncation, truncation_width}`.
- Python (`CEA.pyx`, 5093 LOC): `Reactant`, `Mixture`, `EqSolver/EqSolution/EqDerivatives`,
  `RocketSolver(soln, weights, pc, pi_p, subar, supar, iac=True, n_frz, hc, tc, mdot, ac_at,
  tc_est)`, `ShockSolver(soln, weights, T0, p0, u1|Mach1, reflected, incident_frozen,
  reflected_frozen)`, `DetonationSolver(soln, weights, T1, p1, frozen=False)`, free fns
  `eq_solve/rocket_solve/shock_solve/detonation_solve`, `cea.init`, `cea.units`, `cea.constants.R`.
- Property enums: equilibrium 20 (T,P,V,rho,M,MW,H,U,G,S,gamma_s,frozen/eq Cp,CV,viscosity,
  conductivity, Prandtl); rocket 27 adds `MACH, SONIC_VELOCITY, AE_AT, C_STAR,
  COEFFICIENT_OF_THRUST, ISP, ISP_VACUUM`; shock adds station ratios (`P21/P52, T21/T52,
  RHO12/RHO52, M21/M52` = *molecular-weight* ratios, `V2, U5_P_V2`); detonation 20 (`H1,M1,GAMMA1,
  V_SONIC1, P_P1, T_T1, M_M1`, …). Soft failure `CEA_LAST_VALID_SOLUTION`.
- Other surfaces: legacy CLI `./cea <input>` + `--compile-thermo/--compile-trans`; MATLAB via Python
  bridge; Excel UDF/VBA.
- NASA Python policy: *"No unit conversions are performed by default; inputs and outputs are in the
  documented CEA units."*

## 5. Units + data

- Internal: `R = 8314.51 J/(kmol·K)`; amounts **kmol per kg mixture** (mass basis); P in **bar**,
  T in **K**; properties as `cp/R`, `h/RT`, `s/R` (NASA 9-coefficient convention).
- Output labels (`main.f90`): `P, bar`, `T, K`, `H/U/G, kJ/kg`, `u m/s`, `C* m/s` (+ft/s),
  `Isp m/s` (+lb-s/lb), `Ae/At`.
- Conversions (`units.f90`): k/c/f/r; bar/atm/psi/mmhg; J/kJ/cal/kcal per mole; g/cm³↔kg/m³;
  cm³/g↔m³/kg (unknown units log error and pass through). Python `cea.units` adds m/s↔ft/s.
- Databases: `thermo.inp` 15 802 lines, over 2000 species/phase/reactant records (not all distinct products); `trans.inp` 107 species; compiled
  to `thermo.lib`/`trans.lib`; lookup cwd → `CEA_DATA_DIR` → `./data` → install dir.
- Temperature grids: thermo `200–1000–6000 K` and `1000–6000–20000 K` (`thermo.f90:563-564`);
  transport fits `200–1000`, `1000–5000`, `5000–15000 K`.
- SHA-256: NASA `thermo.inp` `fa774657…84363`, with `\r` stripped `8ad4ca4a…4363` = thermo-species
  copy (identical modulo CRLF); NASA `trans.inp` `379c514a…34d41` = aselimov copy (byte-identical);
  aselimov `thermo.inp` `7a9ada73…fc346` differs from NASA at **line 15782** (one species row) ⇒
  use NASA's file as source of truth.

## 6. Upstream scientific envelope (docs/source/theory.rst, cites RP-1311 Pt1)

- `PV = n_g RT` ⇒ **condensed species occupy zero volume**; upstream warns there are "no guardrails".
- Gibbs minimization with element constraints; Newton + damped update over `[Δln n_j, Δn_j(cond),
  Δln n, Δln T]`; initial guess `n_j = 0.1/NG`, `n = 0.1 kmol/kg`, `T = 3800 K`.
- Trace cutoff `ln n_j < ln(1e-8) = −18.420681` (parameter `trace`); ln kept so species can re-enter.
- Convergence tested on: gas/condensed updates, total moles, T, element balance, entropy, modified
  Lagrange multipliers `π_i`; element residuals scaled by that element's own inventory when
  `> 1e-6 kmol/kg`.
- Condensed activation loop: gas-only solve → add single best G-lowering species → re-converge.
- Derivatives use the activation threshold saved at convergence; optional smooth logistic truncation
  (`truncation_width = 0.25`). Upstream docs leave rocket/shock/detonation/transport/derivatives
  sections as stubs ⇒ validate against RP-1311 sample outputs, not docs.

## 7. Existing Rust work (depth + reuse/license barriers)

**bit-mango/cea-rs — the only real Rust solver (recommended base).**
`MIT OR Apache-2.0` (repo + crates.io `cea-rs 0.1.0`, 2026-08-10). `equilibrium.rs` 3247 LOC +
`mixture.rs` 434 LOC: Gibbs minimization (mass basis, ln n_j + condensed n_j, element reduction for
singular systems) over `nalgebra`. Implements all six modes `EquilibriumMode::{TP,HP,SP,TV,UV,SV}`
and a `Mixture` builder (`with_adiabatic_flame`, `with_heat_addition`,
`with_equilibrium_expansion` = shifting, `with_frozen_expansion` = frozen, `mix`, `scale`).
Controls: `include_condensed`, `include_ions`, `only`, `insert`.
Verified here: `cargo check` clean; `cargo test` **19/19 passed** (1.67 s), test modules keyed to
RP-1311 guide sections (d1, d5, e_hp/e_sp/e_sv/e_tv/e_uv, f_condensed, g_singular, h_ions,
i_postprocess, cea_validation, tp_stability).
Barriers: `thermo-species` declared as **path dep `../thermo-species`** (bare clone fails to build →
switch to registry `0.1` or vendor); no `omit`, no user `trace`, no `smooth_truncation`, no
transport, no CI, `0.1.x` API churn, `panic!` on ions-not-enabled, no H2O(L)/H2O(cr) phase choice.

**thermo-species — required dependency, exact NASA data.**
`MIT OR Apache-2.0`, crates.io `thermo-species 0.1.0` (checksum `1d6b1a27…4216e4`). `build.rs`
parses `data/thermo.inp` (override `THERMO_INP_PATH`) and code-generates the `Species` enum +
`cp/h/s`. Its `thermo.inp` is **byte-identical to NASA's modulo CRLF** (§5). Caveat: build silently
**skips 11 species** with inverted ranges (`Br2(cr)`, `Ca(a)`, `CrN(cr)`, `FeCL3(cr)`, `FeOCL(cr)`,
`Fe3O4(cr)`, `Li(cr)`, `NH4F(cr)`, `Si(cr)`, `Ti3O5(a)`, `U3O8(II)`) — record as known deviation.
No transport data.

**aselimov/cea-rs — data layer only, no solver.** `MIT`. 1316 LOC: `ThermoDB::parse`,
`cp_over_r/h_over_rt/s_over_r`, `TransportDB::parse`, `viscosity_at/conductivity_at`. README
checklist (TP/HP/SP, TV/UV/SV, CJ, shock, rocket) is entirely unchecked and its advertised
`equilibrium/rocket/cj/shock` directories do not exist; `main.rs` is "Hello, world!". Reuse value:
MIT transport parsing (~525 LOC) to close the `trans.inp` gap.

**mcurrie99/nasa-cea — excluded.** FFI wrapper whose `build.rs` CMake/gfortran-builds a *vendored
NASA Fortran tree* (`crates/cea-sys/vendor/cea`); exactly the "Fortran wrapper" excluded by scope.
**No repository license file** (GitHub reports none) ⇒ wrapper code default all-rights-reserved;
do not copy from it.

**Rejected candidates:** crates.io `cea 0.0.0` (Jiri Kovar, MIT) = `println!("Hello, world!")`,
repo `github.com/kovar/cea` returns 404 (tarball kept in `cache/`); `KiThe` (MIT) is a general
kinetics suite, not a CEA port; GitHub search for Rust Gibbs/CEA ports returns only the three above.

## 8. Gap inventory (NASA ⇄ bit-mango)

| Capability | Status |
|---|---|
| TP/HP/SP/TV/UV/SV equilibrium | ✅ implemented + RP-1311-keyed tests |
| Condensed species activation | ✅ active-set and sibling-phase selection; broad coexistence behavior unproven |
| Ions / charge neutrality | ✅ (⚠️ `panic!` instead of `Result`) |
| `only`, `insert` | ✅ ; `omit` ❌ |
| `trace` level, `smooth_truncation`, `truncation_width` | ❌ hardcoded |
| **Rocket `rkt`**: `Ae/At`, `C*`, `Cf`, `Isp`/`Isp_vac`, throat, `iac`, `n_frz`, `mdot`, finite-area, frozen-vs-shifting loop | ❌ only shifting/frozen *expansion states* exist |
| **Shock `shk`** (incident/reflected, frozen variants, stations 1/2/5) | ❌ |
| **Detonation `det`** (CJ, frozen) | ❌ |
| Transport (viscosity/conductivity/Prandtl) | ❌ (aselimov MIT code reusable; needs exact `trans.inp`) |
| Analytic derivatives / sensitivities | ❌ public sensitivities absent; internal cp_eq/cv_eq/gamma_s are thermodynamic responses, not a full sensitivity API |
| O/F ↔ φ ↔ r_eq ↔ fuel% helpers, reactant T-range API | ❌ |
| Legacy `.inp`/`.out`, CLI, `.lib` compile | ❌ (needed only for CLI parity) |
| Public return shape | ⚠️ `{products, T, P, V}` only; NASA solutions expose 20/27/20-property getters + per-station arrays + `converged` |
| Differential validation vs NASA `.out` | ❌ (RP-1311 samples are the oracle) |

## 9. Recommendation for the parent implementer

1. Reuse `bit-mango/cea-rs` + `thermo-species` (only genuine Rust solver; permissive dual license;
   compiles and passes 19 tests here). Fix the `path` dependency to registry/vendored.
2. Vendor NASA `thermo.inp` exactly (SHA §5) with `NOTICE.txt` attribution; never aselimov's copy.
3. Port aselimov's MIT transport parsing if transport is in scope; take NASA `trans.inp` byte-exact.
4. Rocket performance is new code — NASA `rocket.f90` is Apache-2.0 so portable with attribution,
   but its docs section is a stub ⇒ validate against RP-1311 sample outputs.
5. Do not reuse `mcurrie99/nasa-cea` (Fortran wrapper + unlicensed).
6. Envelope: thermo data 200–20 000 K, transport to 15 000 K, zero-volume condensed assumption,
   `initial truncation near 1e-8`, `T_init = 3800 K` (§6); a narrower application envelope must be declared as
   application-specific, not upstream-validated.
7. Ship a THIRD-PARTY/NOTICE file carrying NASA Apache-2.0 + `MIT OR Apache-2.0` (cea-rs,
   thermo-species) + MIT (any aselimov code) with the SHAs from §1.

## 10. Hygiene

No source file modified in any clone; no commits created; `git status --porcelain` clean for
`nasa-cea`, `aselimov-cea-rs`, `bit-mango-cea-rs` (only untracked `._*` AppleDouble files and
gitignored `target/`/`Cargo.lock` from the verification run). All downloads under
`/Volumes/T7 Shield/keplr-work/{upstream,cache}`.
