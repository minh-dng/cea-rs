//! NASA database compilation rules derived from Apache-2.0 NASA CEA database_compile.f90.
// Modified: accept CRLF NASA input; fixed bundled data only. See NOTICE.
use std::collections::{BTreeSet, HashMap};
use std::fmt::Write;
use std::{env, fs, path::Path};

#[path = "src/parsing.rs"]
mod parsing;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rerun-if-changed=data/thermo.inp");
    println!("cargo::rerun-if-changed=src/parsing.rs");
    println!("cargo::rerun-if-changed=src/types.rs");
    let thermo_inp_contents = fs::read_to_string("data/thermo.inp")?;
    let thermo_inp_cleaned: Vec<String> = thermo_inp_contents
        .lines() // Convert file into vec of full lines
        .filter(|row| !row.starts_with("!")) // Remove comments
        .map(|line| line.to_string())
        .collect();

    let mut iter = thermo_inp_cleaned.into_iter().peekable();

    loop {
        if iter
            .peek()
            .ok_or("reached end of input before finding first species (\"e-\")")?
            .split(" ")
            .collect::<Vec<&str>>()[0]
            == "e-"
        {
            break;
        } else {
            iter.next(); // Throw out every line until we hit the first specimen.
        }
    }

    let mut buffer = String::new();
    writeln!(
        &mut buffer,
        "pub use crate::types::{{SpeciesData, AnySpeciesData, ParseSpeciesError}};"
    )?;
    let mut species: HashMap<String, (String, parsing::Metadata, parsing::TemperatureData)> =
        HashMap::new();
    let mut constituents_set: BTreeSet<String> = BTreeSet::new();
    let mut extended_gases = BTreeSet::new();
    let mut previous_phase = 0;
    let mut compiled_phase = 0;
    while iter
        .peek()
        .ok_or("reached end of input before \"END PRODUCTS\"")?
        != "END PRODUCTS"
    {
        let (
            (cleaned_species_identifier, dirty_species_identifier),
            (t_intervals, constituents, phase, mw, h_formation),
            mut cleaned_temperature_data,
        ) = parsing::process_entry(&mut iter);

        // NASA extends two-interval product gases to 20000 K with Cp linear toward its limit at 1e6 K.
        if phase == 0 && cleaned_temperature_data.len() == 2 {
            extended_gases.insert(dirty_species_identifier.clone());
            let (_, t, a) = &cleaned_temperature_data[1];
            let t = *t;
            let atoms: f64 = constituents
                .iter()
                .filter(|(_, e)| e != "E")
                .map(|(n, _)| n)
                .sum();
            let limit = if atoms > 2.1 {
                3.0 * atoms - 1.75
            } else if atoms > 1.9 {
                4.5
            } else {
                2.5
            };
            let cp = a[0] / t.powi(2)
                + a[1] / t
                + a[2]
                + a[3] * t
                + a[4] * t.powi(2)
                + a[5] * t.powi(3)
                + a[6] * t.powi(4);
            let h = -a[0] / t
                + a[1] * t.ln()
                + a[2] * t
                + a[3] * t.powi(2) / 2.0
                + a[4] * t.powi(3) / 3.0
                + a[5] * t.powi(4) / 4.0
                + a[6] * t.powi(5) / 5.0
                + a[7];
            let s = -a[0] / (2.0 * t.powi(2)) - a[1] / t
                + a[2] * t.ln()
                + a[3] * t
                + a[4] * t.powi(2) / 2.0
                + a[5] * t.powi(3) / 3.0
                + a[6] * t.powi(4) / 4.0
                + a[8];
            let slope = (cp - limit) / (t - 1e6);
            let intercept = cp - t * slope;
            cleaned_temperature_data.push((
                6000.0,
                20000.0,
                vec![
                    0.0,
                    0.0,
                    intercept,
                    slope,
                    0.0,
                    0.0,
                    0.0,
                    h - intercept * t - slope * t * t / 2.0,
                    s - intercept * t.ln() - slope * t,
                ],
            ));
        }

        // Save constituents.
        for (_, constituent) in constituents.iter() {
            constituents_set.insert(constituent.clone());
        }

        // NASA numbers each condensed interval, retaining the last number for a named phase.
        // Keep even reversed/zero-width intervals: later intervals can make the species usable.
        if phase == 0 {
            compiled_phase = 0;
        } else {
            for i in 0..cleaned_temperature_data.len() {
                compiled_phase = if phase > previous_phase {
                    compiled_phase + 1
                } else {
                    i as u8 + 1
                };
            }
        }
        previous_phase = phase;

        // Consolidate species.
        if let Some(entry) = species.get_mut(&dirty_species_identifier) {
            // Exact species already exists, check they are the same, and append temperature data together.
            if entry.1.3 != mw || entry.1.4 != h_formation || entry.1.1 != constituents {
                panic!("Species have the exact same name, but different data!");
            }
            // Preserve source order, including overlapping intervals, as NASA's reader does.
            entry.2.extend(cleaned_temperature_data);
            entry.1.2 = compiled_phase;
        } else {
            species.insert(
                dirty_species_identifier,
                (
                    cleaned_species_identifier,
                    (t_intervals, constituents, compiled_phase, mw, h_formation),
                    cleaned_temperature_data,
                ),
            );
        }
    }

    // Add Constituent enum.
    writeln!(
        &mut buffer,
        "#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]"
    )?;
    write!(&mut buffer, "pub enum Constituent {{")?;
    for constituent in constituents_set.iter() {
        write!(&mut buffer, "{},", constituent)?;
    }
    writeln!(&mut buffer, "}}")?;

    let mut all = vec![];
    let mut condensed_phases: HashMap<String, Vec<String>> = HashMap::new();

    let mut consolidated_species = species
        .drain()
        .map(
            // Swap the dirty and clean identifiers
            |(
                dirty_species_identifier,
                (cleaned_species_identifier, metadata, temperature_data),
            )| {
                (
                    cleaned_species_identifier,
                    (dirty_species_identifier, metadata, temperature_data),
                )
            },
        )
        .collect::<Vec<(
            String,
            (
                String,
                (usize, Vec<(f64, String)>, u8, f64, f64),
                Vec<(f64, f64, Vec<f64>)>,
            ),
        )>>();
    consolidated_species.sort_by(|a, b| a.0.cmp(&b.0));
    for (
        cleaned_species_identifier,
        (
            dirty_species_identifier,
            (_, constituents, phase, mw, h_formation),
            cleaned_temperature_data,
        ),
    ) in consolidated_species.iter()
    {
        // Update condensed phase.
        let key = if *phase == 0 {
            // Key is just cleaned identifier.
            cleaned_species_identifier.clone()
        } else {
            let idx_phase_start = cleaned_species_identifier.rfind("__").unwrap_or_else(|| {
                panic!("Failed to find phase in: {}", cleaned_species_identifier)
            });
            let mut identifier_no_phase = cleaned_species_identifier.clone();
            let _ = identifier_no_phase.split_off(idx_phase_start);
            identifier_no_phase
        };
        // Update condensed_phase
        if let Some(entry) = condensed_phases.get_mut(&key) {
            if *phase != 0 {
                entry.push(cleaned_species_identifier.clone());
            }
        } else {
            let mut phases = vec![];
            if *phase != 0 {
                phases.push(cleaned_species_identifier.clone());
            }
            condensed_phases.insert(key, phases);
        }

        // Add each species to the all vec.
        all.push(cleaned_species_identifier.clone());

        let c = constituents.len();
        let t = cleaned_temperature_data.len();
        // Add beginning of static species.
        write!(
            &mut buffer,
            "pub static {}: SpeciesData<{}, {}> = SpeciesData {{",
            cleaned_species_identifier, c, t
        )?;
        // Add symbol.
        write!(&mut buffer, "symbol: \"{}\",", dirty_species_identifier)?;
        // Add constituents.
        let mut constituent_str = "".to_string();
        for (moles, constituent) in constituents {
            constituent_str += &format!("({:?}, Constituent::{}),", moles, constituent);
        }
        write!(&mut buffer, "constituents: [{}],", constituent_str)?;
        // Add temperature data.
        let mut temperature_data_str = "".to_string();
        for (low_temperature, high_temperature, coefficients) in cleaned_temperature_data {
            let mut coefficients_str = "".to_string();
            for coefficient in coefficients {
                coefficients_str += &format!("{:?},", coefficient);
            }
            temperature_data_str += &format!(
                "({:?}, {:?}, [{}],),",
                low_temperature, high_temperature, coefficients_str
            );
        }
        write!(&mut buffer, "temperature_data: [{}],", temperature_data_str)?;
        // Add mw.
        write!(&mut buffer, "mw: {:?},", mw)?;
        // Add h_formation.
        write!(&mut buffer, "h_formation: {:?},", h_formation)?;
        // Add phase.
        write!(&mut buffer, "phase: {:?},", phase)?;
        let extrapolated_above = if extended_gases.contains(dirty_species_identifier) {
            Some(6000.0)
        } else {
            None
        };
        write!(&mut buffer, "extrapolated_above: {:?},", extrapolated_above)?;
        // Close out struct definition.
        writeln!(&mut buffer, "}};")?;
    }

    // Add Species enum.
    writeln!(
        &mut buffer,
        "#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]"
    )?;
    write!(&mut buffer, "pub enum Species {{")?;
    for species in all.iter() {
        write!(&mut buffer, "{},", species)?;
    }
    writeln!(&mut buffer, "}}")?;

    write!(&mut buffer, "pub static ALL: [Species; {}] = [", all.len())?;
    for species in all.iter() {
        write!(&mut buffer, "Species::{},", species)?;
    }
    writeln!(&mut buffer, "];")?;

    // Start Species impl
    write!(&mut buffer, "impl Species {{")?;

    // Add all
    writeln!(&mut buffer, "pub fn all() -> &'static [Species] {{&ALL}}")?;

    // Add phases
    write!(&mut buffer, "pub fn phases(&self) -> &'static [Species] {{")?;
    write!(&mut buffer, "match self {{")?;
    for (base, phases) in condensed_phases.iter() {
        let base_is_real = all.contains(base);
        let mut left = if base_is_real {
            format!("Species::{}", base)
        } else {
            String::new()
        };
        let mut right = "&[".to_string();
        phases.iter().for_each(|phase| {
            left += &format!(" | Species::{}", phase);
            right += &format!("Species::{},", phase);
        });
        right += "],";
        write!(&mut buffer, "{} => {}", left, right)?;
    }
    // Close match arm
    write!(&mut buffer, "}}")?;
    // Close fn definition
    writeln!(&mut buffer, "}}")?;

    // Add data
    write!(
        &mut buffer,
        "pub fn data(&self) -> &'static dyn AnySpeciesData {{"
    )?;
    write!(&mut buffer, "match self {{")?;
    for species in all.iter() {
        write!(&mut buffer, "Species::{} => &{},", species, species)?;
    }
    // Close match arm
    write!(&mut buffer, "}}")?;
    // Close fn definition
    write!(&mut buffer, "}}")?;

    // Close the implementation.
    writeln!(&mut buffer, "}}")?;

    // Implement from_str for Species.
    writeln!(&mut buffer, "impl std::str::FromStr for Species {{")?;
    writeln!(&mut buffer, "type Err = ParseSpeciesError;")?;
    writeln!(
        &mut buffer,
        "fn from_str(s: &str) -> Result<Self, Self::Err> {{"
    )?;
    write!(&mut buffer, "match s {{")?;
    for (cleaned_species_identifier, (dirty_species_identifier, (_, _, _, _, _), _)) in
        consolidated_species.iter()
    {
        write!(
            &mut buffer,
            "\"{}\" => Ok(Species::{}),",
            dirty_species_identifier, cleaned_species_identifier
        )?;
    }
    write!(&mut buffer, "_ => Err(ParseSpeciesError), }}",)?;

    // Close functin definition.
    writeln!(&mut buffer, "}}")?;

    // Close the implementation
    writeln!(&mut buffer, "}}")?;

    // Define species!() macro.
    writeln!(&mut buffer, "#[macro_export]")?;
    writeln!(&mut buffer, "macro_rules! species {{")?;
    for (cleaned_species_identifier, (dirty_species_identifier, (_, _, _, _, _), _)) in
        consolidated_species.iter()
    {
        write!(
            &mut buffer,
            "(\"{}\") => {{ $crate::Species::{}}};",
            dirty_species_identifier, cleaned_species_identifier
        )?;
    }
    write!(
        &mut buffer,
        "($other:literal) => {{ compile_error!(concat!(\"Unknown species: \", $other))}}",
    )?;
    // Close the macro
    writeln!(&mut buffer, "}}")?;

    let out_dir =
        env::var_os("OUT_DIR").ok_or("OUT_DIR not set — build scripts must be run via cargo")?;
    let dest_path = Path::new(&out_dir).join("thermo.rs");
    fs::write(&dest_path, buffer)?;

    Ok(())
}
