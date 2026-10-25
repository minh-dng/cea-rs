# Source and data provenance

This independent Keplr implementation contains genuine Rust thermochemistry
and equilibrium code. It neither links to Fortran CEA nor calls a CEA process
or Python package at runtime. NASA CEA is used only by the comparison tools.

- NASA CEA: [`2f79a647e85737742fb08a127da00d1ea1d069da`](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da), current main inspected 2026-10-25, package version 3.3.4. Apache-2.0. The original LICENSE and NOTICE are retained at repository root. NASA copyright: 2025 United States Government as represented by the Administrator of NASA.
- [`data/thermo.inp`](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp) is copied without content or line-ending changes into `vendor/thermo-species/data/thermo.inp`. SHA-256: `fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`. `.gitattributes` preserves its bytes. The data are embedded in compiled Rust code; no network download or Fortran-generated binary database is needed after installation. Reference: [NASA TP-2002-211556](https://ntrs.nasa.gov/citations/20020085330).
- Rust equilibrium algorithm: [`bit-mango/cea-rs` at `914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad`](https://github.com/bit-mango/cea-rs/tree/914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad), Copyright (c) 2026 Bit Mango, originally MIT OR Apache-2.0. Used under Apache-2.0; both original license texts retained in `vendor/`. Adapted file: `src/equilibrium_core.rs`.
- Thermodynamic polynomial types and generator: [`bit-mango/thermo-species` at `d4061dddda7f04e2edf75cdb235d19d1fed4161b`](https://github.com/bit-mango/thermo-species/tree/d4061dddda7f04e2edf75cdb235d19d1fed4161b), Copyright (c) 2026 Bit Mango, MIT OR Apache-2.0. Used under Apache-2.0; original licenses retained within the vendored crate. Its source database has the same logical lines as NASA's file, but LF line endings; we replaced that copy with the byte-identical NASA file.
- [`aselimov/cea-rs`](https://github.com/aselimov/cea-rs) was inspected before implementation. It provides thermodynamic/transport parsing but its README marks the equilibrium/application solvers unimplemented. No code was copied from it.

Changes to reused Rust code are marked in the files and NOTICE: arbitrary
reactant formulas; strict input and output checks; fail-closed equilibrium
derivatives; post-update constant-volume entropy closure; fixed bundled data;
CRLF support; NASA's gas constant; Python bindings; IAC nozzle expansion.
The Bit Mango high-level mixture wrapper is not used: it has a different
mole-based interface and its own rounded gas constant.

The generator omits 11 NASA condensed records with empty or reversed declared
temperature intervals: Br2(cr), Ca(a), CrN(cr), FeCL3(cr), FeOCL(cr), Fe3O4(cr),
Li(cr), NH4F(cr), Si(cr), Ti3O5(a), U3O8(II). Other valid phases remain available.
This omission is a limitation, not a repair to NASA's database.

NASA source/API inspection covers `source/bind/python/cea/lib/libcea.pyi`,
`source/fits.f90`, `source/thermo.f90`, `source/equilibrium.f90`,
`source/rocket.f90`, `source/shock.f90`, `source/detonation.f90` and transport
modules. See the implementation matrix for the substantial API and mode gaps.

Consumer seams were read without edits at Hybrid `b249d5d` and its
`feat-liquid-design-point` checkout: `RocketCEA`, custom ABS/paraffin/HTPB,
and liquid thermochemistry. Their production dependency stays native `cea`.
