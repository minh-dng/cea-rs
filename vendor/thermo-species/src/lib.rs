//! A lightweight, type-safe wrapper around NASA's 9-coefficient
//! thermodynamic polynomials — the same species data format NASA CEA
//! (Chemical Equilibrium with Applications) reads from `thermo.inp`.
//!
//! At build time, this crate parses `data/thermo.inp` and code-generates a
//! [`Species`] enum plus one const [`SpeciesData`] per species, so species
//! lookups and property evaluation are checked at compile time wherever
//! possible instead of failing at runtime on a typo'd symbol.
//!
//! **Status: under active development.** This is a `0.1.x` crate — the API
//! may still change between minor versions while it stabilizes.
//!
//! # Getting a `Species`
//!
//! There are two ways, depending on whether the symbol is known at the call
//! site or only discovered at runtime.
//!
//! ## Compile time, via [`species!`]
//!
//! Use this whenever the species you want is written directly in your code.
//! An unknown symbol is a compile error, not a panic or a `Result` to
//! unwrap:
//!
//! ```
//! use thermo_species::species;
//!
//! let water = species!("H2O");
//! let solid_water = species!("H2O(cr)"); // condensed phase
//! ```
//!
//! ## Runtime, via [`Species::from_str`]
//!
//! Use this when the species symbol comes from user input, a config file, or
//! anything else not known until runtime. It returns a `Result` instead of
//! failing to compile:
//!
//! ```
//! use std::str::FromStr;
//! use thermo_species::Species;
//!
//! let symbol = "H2O"; // e.g. read from a file or user input
//! let water = Species::from_str(symbol)?;
//! # Ok::<(), thermo_species::ParseSpeciesError>(())
//! ```
//!
//! # Reading thermodynamic properties
//!
//! Every [`Species`] exposes its data via [`Species::data`], which returns
//! `&dyn `[`AnySpeciesData`]:
//!
//! ```
//! use thermo_species::species;
//!
//! let water = species!("H2O");
//! let data = water.data();
//!
//! let cp = data.cp(1000.0); // J/(mol·K) at 1000 K
//! let h = data.h(1000.0);   // J/mol at 1000 K
//! let s = data.s(1000.0);   // J/(mol·K) at 1000 K
//! let mw = data.mw();       // g/mol
//! ```
//!
//! R-normalized variants ([`AnySpeciesData::cp_over_r`],
//! [`AnySpeciesData::h_over_r`], [`AnySpeciesData::s_over_r`],
//! [`AnySpeciesData::u_over_r`]) are also available — these match NASA CEA's
//! own internal units (Kelvin-scaled, gas-constant-free), which is
//! convenient when porting or cross-checking against CEA output directly.

mod parsing;
mod types;
#[allow(non_upper_case_globals, non_camel_case_types)]
mod generated {
    include!(concat!(env!("OUT_DIR"), "/thermo.rs"));
}
pub use generated::{AnySpeciesData, Constituent, ParseSpeciesError, Species, SpeciesData};

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_species_macro() {
        // Compile time Species
        assert_eq!(species!("H2O(cr)"), Species::H2O__1);
    }
    #[test]
    fn test_from_str() {
        // Run time dynamic Species
        assert_eq!(Species::from_str("H2O").unwrap(), Species::H2O);

        println!("{}", species!("CH4"));
    }
}
