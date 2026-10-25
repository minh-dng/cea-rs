# thermo-species

A lightweight, type-safe wrapper around NASA's 9-coefficient thermodynamic
polynomials (the same species data format used by NASA CEA — Chemical
Equilibrium with Applications). At build time, `thermo-species` parses
`data/thermo.inp` and code-generates a `Species` enum, one const struct per
species, and Cp/H/S evaluation from the polynomial fits — so species lookups
and thermodynamic property calculations are checked at compile time wherever
possible, instead of failing at runtime on a typo'd species name.

> **Status: under active development.** This is a `0.1.x` crate — the API may
> still change between minor versions while it stabilizes. It's published now
> because it's useful as-is, not because it's finished.

## Two ways to get a `Species`

**1. Compile time, via the `species!` macro** — use this whenever the species
you want is known at the call site. An unknown symbol is a compile error, not
a runtime one:

```rust
use thermo_species::species;

let water = species!("H2O");
let carbon_dioxide = species!("CO2");
let solid_water = species!("H2O(cr)"); // condensed phase

// species!("H2Oo") would fail to compile:
// error: Unknown species: "H2Oo"
```

**2. Runtime, via `Species::from_str`** — use this when the species symbol
comes from user input, a config file, or anywhere else not known until
runtime. It returns a `Result` instead of failing to compile:

```rust
use std::str::FromStr;
use thermo_species::Species;

let symbol = "H2O"; // e.g. read from a file or user input
let water = Species::from_str(symbol)?;
# Ok::<(), thermo_species::ParseSpeciesError>(())
```

## Reading thermodynamic properties

Every `Species` exposes its data via `.data()`, which returns
`&dyn AnySpeciesData`:

```rust
use thermo_species::species;

let water = species!("H2O");
let data = water.data();

let cp = data.cp(1000.0);   // J/(mol·K) at 1000 K
let h = data.h(1000.0);     // J/mol at 1000 K
let s = data.s(1000.0);     // J/(mol·K) at 1000 K
let mw = data.mw();         // g/mol

println!("{water}"); // Display impl: symbol, phase, mw, ΔHf, constituents, ranges
```

R-normalized variants (`cp_over_r`, `h_over_r`, `s_over_r`, `u_over_r`) are
also available — these match NASA CEA's own internal units (Kelvin-scaled,
gas-constant-free), which is convenient when porting or cross-checking
against CEA output directly.

## Condensed phases

Species with multiple phases (e.g. `H2O` gas vs. `H2O(cr)` solid vs. `H2O(L)`
liquid) are distinct `Species` variants. `.phases()` lists the sibling phases
of a given species:

```rust
use thermo_species::species;

let water_gas = species!("H2O");
for phase in water_gas.phases() {
    println!("{phase}");
}
```

## Data source

`data/thermo.inp` is NASA's published 9-coefficient polynomial dataset (the
same input file CEA itself reads). A `trans.inp`-derived transport-properties
database (viscosity/thermal conductivity coefficients) is a possible future
addition, but isn't included yet.

## License

Licensed under either of

- MIT license ([LICENSE-MIT](LICENSE-MIT))
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))

at your option.
