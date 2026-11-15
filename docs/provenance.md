# Source and data provenance

This independent Keplr implementation contains genuine Rust thermochemistry
and equilibrium code. It neither links to Fortran CEA nor calls a CEA process
or Python package at runtime. NASA CEA is used only by the comparison tools.

- NASA CEA: [`2f79a647e85737742fb08a127da00d1ea1d069da`](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da), current main inspected 2026-10-25, package version 3.3.4. Apache-2.0. The original LICENSE and NOTICE are retained at repository root. NASA copyright: 2025 United States Government as represented by the Administrator of NASA.
- `src/detonation.rs` derives its CJ initialization, damped Newton iteration and stopping test from NASA's `source/detonation.f90` at that commit (Apache-2.0). Volume derivatives use public equilibrium solves and the Cp/gamma identity. Returned states are re-evaluated and independently checked for jump conservation and sonic closure; native transport calls and unfinished frozen detonation are excluded.
- [`data/thermo.inp`](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp) is copied without content or line-ending changes into `vendor/thermo-species/data/thermo.inp`. SHA-256: `fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`. `.gitattributes` preserves its bytes. The data are embedded in compiled Rust code; no network download or Fortran-generated binary database is needed after installation. Reference: [NASA TP-2002-211556](https://ntrs.nasa.gov/citations/20020085330).
- [`data/trans.inp`](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/trans.inp) is copied byte-for-byte into `vendor/nasa-transport/trans.inp`. SHA-256: `379c514a7f7638371d8d9254a1d653cd078084a9371061aa0b45c5a5c3334d41`. `.gitattributes` preserves its CRLF and padding. The source distribution includes it and wheels embed it with `include_str!`; the root Apache-2.0 LICENSE and NOTICE cover this NASA data, as for `thermo.inp`.
- `src/transport.rs` is an Apache-2.0 derived work of NASA `source/transport.f90`, `source/database_compile.f90` (text format), and the transport basis, mixing and reaction calculations in `source/equilibrium.f90` at the same pinned commit. It uses NASA's rounded constants and end-fit extrapolation. Native `cea==3.3.4`, initialized with `trans=True`, supplied the 40 goldens and generation instructions in `tests/transport.rs` on 2026-11-15. All five SI outputs pass a relative tolerance of `1e-4`; the worst measured difference is `2.81629e-5` (CH4/Air at 2000 K, equilibrium conductivity). Cases cover H2/O2 and CH4/Air at 500–5000 K, all three RP-1311 example 2 TV states, all four example 6 burned states, the example 11 chamber TP state with ions, graphite, ionized hydrogen, argon interval boundaries/extrapolation, and frozen-at-chamber/throat nozzle states. As in NASA, the equilibrium transport columns on a frozen composition are reaction-model diagnostics. Example 11's native chamber transport matches at the native TP state to `2.18e-8` relative; the example 11 rocket stations themselves match native (see `docs/rocket-parity/README.md`).
- Rust equilibrium algorithm: [`bit-mango/cea-rs` at `914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad`](https://github.com/bit-mango/cea-rs/tree/914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad), Copyright (c) 2026 Bit Mango, originally MIT OR Apache-2.0. Used under Apache-2.0; both original license texts retained in `vendor/`. Adapted file: `src/equilibrium_core.rs`.
- Thermodynamic polynomial types and generator: [`bit-mango/thermo-species` at `d4061dddda7f04e2edf75cdb235d19d1fed4161b`](https://github.com/bit-mango/thermo-species/tree/d4061dddda7f04e2edf75cdb235d19d1fed4161b), Copyright (c) 2026 Bit Mango, MIT OR Apache-2.0. Used under Apache-2.0; original licenses retained within the vendored crate. Its source database has the same logical lines as NASA's file, but LF line endings; we replaced that copy with the byte-identical NASA file.
- [`aselimov/cea-rs`](https://github.com/aselimov/cea-rs) was inspected before implementation. It provides thermodynamic/transport parsing but its README marks the equilibrium/application solvers unimplemented. No code was copied from it.
- Mixture conversions, atomic valences, temperature-range queries and equilibrium sensitivities in `src/mixture.rs`, `src/thermo.rs` and `src/derivatives.rs` are derived from NASA CEA `source/mixture.f90`, `source/atomic_data.f90`, `source/thermo.f90` and `source/equilibrium.f90` at the NASA commit above, under Apache-2.0. Native `cea==3.3.4` golden generation and SI conversions are recorded in `tests/mixture.rs` and `tests/derivatives.rs`.
- Shock equations, initialization, damping and reflected-frozen reporting in `src/shock.rs` are derived from NASA [`source/shock.f90`](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/shock.f90), under Apache-2.0. The Rust adaptation uses public equilibrium solves with numerical TP derivatives, excludes transport and rejects retained/nonconverged states. `tests/shock.rs` records native `cea==3.3.4` goldens and their generation procedure.

Changes to reused Rust code are marked in the files and NOTICE: arbitrary
reactant formulas; strict input and output checks; fail-closed equilibrium
derivatives; post-update constant-volume entropy closure; fixed bundled data;
CRLF support; NASA's gas constant; Python bindings; IAC nozzle expansion.
The Bit Mango high-level mixture wrapper is not used: it has a different
mole-based interface and its own rounded gas constant.

The 11 NASA condensed records with empty or reversed declared temperature
intervals (Br2(cr), Ca(a), CrN(cr), FeCL3(cr), FeOCL(cr), Fe3O4(cr), Li(cr),
NH4F(cr), Si(cr), Ti3O5(a), U3O8(II)) are included.

NASA source/API inspection covers `source/bind/python/cea/lib/libcea.pyi`,
`source/fits.f90`, `source/thermo.f90`, `source/equilibrium.f90`,
`source/rocket.f90`, `source/shock.f90`, `source/detonation.f90` and transport
modules. See the implementation matrix for the substantial API and mode gaps.

Consumer seams were read without edits at Hybrid `b249d5d` and its
`feat-liquid-design-point` checkout: `RocketCEA`, custom ABS/paraffin/HTPB,
and liquid thermochemistry. Their production dependency stays native `cea`.

- Equilibrium parity changes derive from Apache-2.0 NASA CEA at the commit above: `equilibrium.f90`, `thermo.f90`, `mixture.f90` and `database_compile.f90`. Adapted modules: `src/equilibrium_core.rs`, `src/thermo.rs` and `vendor/thermo-species/{build.rs,src/types.rs}`. This supersedes the older omission note above: all 11 condensed records are now retained, with NASA interval order and compiled phase numbering. Source database bytes and NOTICE are unchanged.
- `src/rocket.rs`: Apache-2.0 derived work ported from NASA CEA `source/rocket.f90` at `2f79a647e85737742fb08a127da00d1ea1d069da`: IAC/FAC stations, pressure/area Newton iterations, station-index freezing, frozen fit limits and ideal performance; original NASA NOTICE remains intact.
