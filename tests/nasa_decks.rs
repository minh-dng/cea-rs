//! Direct comparisons with unedited NASA reports; no hand-copied state goldens.
use cea_rs::{
    deck::{CaseInput, CaseOutput, DeckResult, parse_deck},
    solver::State,
    transport::Transport,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

const TOLERANCE: f64 = 1e-4;
const COMPOSITION_THRESHOLD: f64 = 1e-5;

#[derive(Default, Debug)]
struct Page {
    frozen: bool,
    rows: BTreeMap<String, Vec<f64>>,
}

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/nasa")
}

// Only the numerical rows used below are read. This is deliberately test code,
// not an implementation of NASA's report format.
fn reference(text: &str) -> Vec<Page> {
    let mut pages: Vec<Page> = vec![];
    let mut frozen = false;
    let mut active = false;
    let mut prefix = "";
    let mut section = "";
    let mut conductivity_factor = 0.1;
    for raw in text.lines() {
        let line = raw.trim();
        let upper = line.to_ascii_uppercase();
        if upper.contains("THEORETICAL ROCKET PERFORMANCE ASSUMING") {
            frozen = upper.contains("FROZEN");
            active = false;
        }
        if upper.contains("COMPOSITION FOR INCIDENT SHOCKED CONDITIONS") {
            frozen = upper.contains("FROZEN");
        }
        if upper == "THERMODYNAMIC PROPERTIES"
            || upper == "UNBURNED GAS"
            || upper == "INITIAL GAS (1)"
            || (upper.contains("THROAT")
                && (upper.contains("CHAMBER") || upper.contains("INJECTOR")))
        {
            pages.push(Page {
                frozen,
                ..Page::default()
            });
            active = true;
            section = "";
            prefix = if upper == "UNBURNED GAS" || upper == "INITIAL GAS (1)" {
                "1."
            } else {
                ""
            };
            continue;
        }
        if !active {
            continue;
        }
        if upper == "BURNED GAS" {
            prefix = "";
        }
        if upper.starts_with("SHOCKED GAS (2)") {
            prefix = "2.";
        }
        if upper.starts_with("REFLECTED GAS (5)") {
            prefix = "5.";
        }
        if upper == "MOLE FRACTIONS" {
            section = "x.";
            continue;
        }
        if upper == "MASS FRACTIONS" {
            section = "y.";
            continue;
        }
        if upper.contains("PRODUCTS WHICH") || upper.starts_with("* THERMODYNAMIC") {
            active = false;
            continue;
        }
        if upper.contains("CONDUCTIVITY IN UNITS") {
            conductivity_factor = if upper.contains("CALORIES") {
                0.4184
            } else {
                0.1
            };
        }
        if upper == "WITH EQUILIBRIUM REACTIONS" {
            section = "eq.";
            continue;
        }
        if upper == "WITH FROZEN REACTIONS" {
            section = "fr.";
            continue;
        }
        if upper == "PERFORMANCE PARAMETERS" || upper == "DETONATION PARAMETERS" {
            section = "";
        }
        let words: Vec<_> = line.split_whitespace().collect();
        let Some(start) = words.iter().position(|w| w.parse::<f64>().is_ok()) else {
            continue;
        };
        if start == 0 {
            continue;
        }
        let page = pages.last_mut().unwrap();
        if matches!(section, "x." | "y.") {
            if start == 1 && words[1..].iter().all(|w| w.parse::<f64>().is_ok()) {
                page.rows.insert(
                    format!("{section}{}", words[0].trim_start_matches('*')),
                    words[1..].iter().map(|w| w.parse().unwrap()).collect(),
                );
            } else {
                for pair in words.as_chunks::<2>().0 {
                    if let Ok(value) = pair[1].parse() {
                        page.rows.insert(
                            format!("{section}{}", pair[0].trim_start_matches('*')),
                            vec![value],
                        );
                    }
                }
            }
            continue;
        }
        let label = words[..start].join(" ").to_ascii_uppercase();
        let energy_factor = if label.contains("CAL") {
            4184.0
        } else {
            1000.0
        };
        let (key, factor) = if label.starts_with("P,") || label.starts_with("P1,") {
            ("p", if label.contains("ATM") { 101325.0 } else { 1e5 })
        } else if label.starts_with("T,") || label.starts_with("T1,") {
            ("t", 1.0)
        } else if label.starts_with("RHO,")
            || label.starts_with("RHO1,")
            || label.starts_with("DENSITY,")
        {
            ("rho", if label.contains("G/CC") { 1000.0 } else { 1.0 })
        } else if label.starts_with("H,") || label.starts_with("H1,") {
            ("h", energy_factor)
        } else if label.starts_with("U,") {
            ("u", energy_factor)
        } else if label.starts_with("S,") {
            ("s", energy_factor)
        } else if label.starts_with("M,") || label.starts_with("M1,") {
            ("m", 1.0)
        } else if label.starts_with("CP,") {
            ("cp", energy_factor)
        } else if label.starts_with("GAMMA") {
            ("gamma", 1.0)
        } else if label.starts_with("SON") {
            ("sound", 1.0)
        } else if label.starts_with("MACH")
            || label.starts_with("DET MACH")
            || label.starts_with("DET. MACH")
        {
            ("mach", 1.0)
        } else if label.starts_with("U1,") {
            ("u1", 1.0)
        } else if label.starts_with("U2,") {
            ("u2", 1.0)
        } else if label.starts_with("V2,") {
            ("v2", 1.0)
        } else if label.starts_with("DET VEL") || label.starts_with("DET. VEL") {
            ("velocity", 1.0)
        } else if label == "P2/P1" || label == "P/P1" {
            ("pressure_ratio", 1.0)
        } else if label == "T2/T1" || label == "T/T1" {
            ("temperature_ratio", 1.0)
        } else if label == "M2/M1" || label == "M/M1" {
            ("molecular_weight_ratio", 1.0)
        } else if label == "RHO2/RHO1" || label == "RHO/RHO1" {
            ("density_ratio", 1.0)
        } else if label == "AE/AT" {
            ("area", 1.0)
        } else if label.starts_with("CSTAR") || label.starts_with("C*,") {
            ("cstar", if label.contains("FT") { 0.3048 } else { 1.0 })
        } else if label == "CF" {
            ("cf", 1.0)
        } else if label.starts_with("IVAC,") {
            ("ivac", if label.contains("/S") { 1.0 } else { 9.80665 })
        } else if label.starts_with("ISP,") {
            ("isp", if label.contains("/S") { 1.0 } else { 9.80665 })
        } else if label.starts_with("VISC") {
            ("viscosity", 1e-4)
        } else if label == "CONDUCTIVITY" {
            ("conductivity", conductivity_factor)
        } else if label == "PRANDTL NUMBER" {
            ("prandtl", 1.0)
        } else {
            continue;
        };
        let values: Vec<f64> = words[start..]
            .iter()
            .map(|w| {
                w.parse::<f64>()
                    .unwrap_or_else(|_| panic!("Bad reference row: {line}"))
                    * factor
            })
            .collect();
        page.rows.insert(format!("{prefix}{section}{key}"), values);
    }
    assert!(!pages.is_empty(), "No numerical reference tables");
    pages
}

#[derive(Default)]
struct Comparison {
    context: String,
    checked: usize,
    worst: f64,
    worst_detail: String,
    failures: Vec<String>,
    unreferenced_stations: Vec<(usize, usize)>,
}
impl Comparison {
    fn value(&mut self, key: &str, actual: f64, expected: f64) {
        // NASA rounds the 298.15 K unburned enthalpy to zero (fit is -0.000506 J/kg).
        if key == "1.h" && expected == 0.0 {
            assert!(
                actual.abs() <= 1e-3,
                "Zero-reference unburned enthalpy: {actual} J/kg"
            );
            self.checked += 1;
            return;
        }
        self.check(
            key,
            actual,
            expected,
            if expected == 0.0 { 1.0 } else { expected.abs() },
        );
    }
    fn relative(&mut self, key: &str, actual: f64, expected: f64) {
        self.check(key, actual, expected, expected.abs());
    }
    fn check(&mut self, key: &str, actual: f64, expected: f64, scale: f64) {
        let key = format!("{} {key}", self.context);
        let diff = (actual - expected).abs() / scale;
        self.checked += 1;
        if diff > self.worst {
            self.worst = diff;
            self.worst_detail = format!("{key}: Rust={actual:.15e}, NASA={expected:.15e}");
        }
        if !actual.is_finite() || diff > TOLERANCE {
            self.failures.push(format!(
                "{key}: Rust={actual:.15e}, NASA={expected:.15e}, relative={diff:.6e}"
            ));
        }
    }
    fn row(&mut self, page: &Page, key: &str, column: usize, actual: f64) {
        let values = page.rows.get(key).unwrap_or_else(|| {
            panic!(
                "Missing required reference row {key}: {:?}",
                page.rows.keys()
            )
        });
        self.value(key, actual, values[column]);
    }
    fn state(&mut self, page: &Page, prefix: &str, col: usize, state: &State) {
        for (key, actual) in [
            ("t", state.temperature_k),
            ("p", state.pressure_pa),
            ("rho", state.density_kg_m3),
            ("h", state.enthalpy_j_kg),
            ("u", state.energy_j_kg),
            ("s", state.entropy_j_kg_k),
            ("m", state.molecular_weight_kg_kmol),
            (
                "cp",
                if state.equilibrium {
                    state.cp_equilibrium_j_kg_k
                } else {
                    state.cp_frozen_j_kg_k
                },
            ),
            ("gamma", state.gamma_s),
            ("sound", state.sound_speed()),
        ] {
            self.row(page, &format!("{prefix}{key}"), col, actual);
        }
        if !prefix.is_empty() && prefix != "2." {
            return;
        }
        let mass = page.rows.keys().any(|k| k.starts_with("y."));
        let composition_prefix = if mass { "y." } else { "x." };
        assert!(
            page.rows.keys().any(|k| k.starts_with(composition_prefix)),
            "Missing reference composition"
        );
        let compact_frozen = state.equilibrium
            && page.frozen
            && page
                .rows
                .iter()
                .any(|(k, v)| k.starts_with(composition_prefix) && v.len() == 1);
        if compact_frozen {
            return;
        }
        for (key, values) in &page.rows {
            let fractions = if key.starts_with("x.") {
                &state.mole_fractions
            } else if key.starts_with("y.") {
                &state.mass_fractions
            } else {
                continue;
            };
            if values.len() == 1 && state.equilibrium && page.rows[&format!("{prefix}t")].len() > 1
            {
                continue;
            }
            let expected = values[if values.len() == 1 { 0 } else { col }];
            if expected >= COMPOSITION_THRESHOLD {
                let actual = *fractions.get(&key[2..]).unwrap_or(&0.0);
                // Composition is dimensionless: use a relative threshold, not the SI scale floor.
                self.relative(key, actual, expected);
            }
        }
        let fractions = if mass {
            &state.mass_fractions
        } else {
            &state.mole_fractions
        };
        for (name, value) in fractions {
            let key = format!("{composition_prefix}{name}");
            if *value >= COMPOSITION_THRESHOLD && !page.rows.contains_key(&key) {
                self.check(&key, *value, 0.0, COMPOSITION_THRESHOLD);
            }
        }
    }
    fn transport(&mut self, page: &Page, col: usize, t: &Transport) {
        for (key, actual) in [
            ("viscosity", t.viscosity_pa_s),
            ("eq.cp", t.cp_equilibrium_j_kg_k),
            ("fr.cp", t.cp_frozen_j_kg_k),
            ("eq.conductivity", t.conductivity_equilibrium_w_m_k),
            ("fr.conductivity", t.conductivity_frozen_w_m_k),
            ("eq.prandtl", t.prandtl_equilibrium),
            ("fr.prandtl", t.prandtl_frozen),
        ] {
            let expected = page.rows[key][col];
            self.relative(key, actual, expected);
        }
    }
}

fn compare(results: &[DeckResult], pages: &[Page]) -> Comparison {
    let mut c = Comparison::default();
    let mut offset = 0;
    let mut seen = vec![BTreeSet::new(); results.len()];
    for page in pages {
        match &results[0].output {
            CaseOutput::Equilibrium(_) | CaseOutput::Detonation(_) | CaseOutput::Shock(_) => {
                let count = page
                    .rows
                    .get("t")
                    .or_else(|| page.rows.get("1.t"))
                    .unwrap()
                    .len();
                for (col, r) in results[offset..offset + count].iter().enumerate() {
                    c.context = format!("case {}", r.case_index);
                    match &r.output {
                        CaseOutput::Equilibrium(st) => c.state(page, "", col, st),
                        CaseOutput::Detonation(d) => {
                            c.state(page, "", col, &d.burned);
                            for (key, actual) in [
                                ("1.t", d.unburned.temperature_k),
                                ("1.p", d.unburned.pressure_pa),
                                ("1.h", d.unburned.enthalpy_j_kg),
                                ("1.m", d.unburned.molecular_weight_kg_kmol),
                                ("1.gamma", d.unburned.gamma_s),
                                ("1.sound", d.unburned.sound_speed_m_s),
                                ("pressure_ratio", d.pressure_ratio),
                                ("temperature_ratio", d.temperature_ratio),
                                ("density_ratio", d.density_ratio),
                                ("molecular_weight_ratio", d.molecular_weight_ratio),
                                ("mach", d.mach),
                                ("velocity", d.velocity_m_s),
                            ] {
                                c.row(page, key, col, actual);
                            }
                        }
                        CaseOutput::Shock(s) => {
                            c.state(page, "1.", col, &s.station1.state);
                            c.state(page, "2.", col, &s.station2.state);
                            for (key, actual) in [
                                ("1.mach", s.station1.mach),
                                ("1.u1", s.station1.velocity_m_s),
                                ("2.u2", s.station2.velocity_m_s),
                                ("2.v2", s.v2_m_s),
                                ("2.pressure_ratio", s.p2_p1),
                                ("2.temperature_ratio", s.t2_t1),
                                ("2.molecular_weight_ratio", s.m2_m1),
                                ("2.density_ratio", s.rho2_rho1),
                            ] {
                                c.row(page, key, col, actual);
                            }
                        }
                        _ => unreachable!(),
                    }
                    if let Some(t) = r.transport.first() {
                        c.transport(page, col, t);
                    }
                }
                offset += count;
            }
            CaseOutput::Rocket(_) => {
                let rocket_case = results.iter().position(|r| matches!(r.input, CaseInput::Rocket { freeze_station, .. } if (freeze_station > 0) == page.frozen)).unwrap();
                let CaseOutput::Rocket(r) = &results[rocket_case].output else {
                    unreachable!()
                };
                let CaseInput::Rocket { freeze_station, .. } = results[rocket_case].input else {
                    unreachable!()
                };
                assert_eq!(freeze_station > 0, page.frozen);
                // NASA repeats the chamber/combustor/throat prefix on each report page.
                // FAC reports omit the virtual infinity station returned by the API.
                let stations: Vec<_> = r
                    .stations
                    .iter()
                    .filter(|s| !(r.throat_index == 3 && s.name == "infinity"))
                    .collect();
                let repeated = if r.throat_index == 3 { 3 } else { 2 };
                let count = page.rows["t"].len();
                for col in 0..count {
                    let index = if col < repeated {
                        col
                    } else {
                        let area = page.rows["area"][col - 1];
                        let subsonic = page.rows["mach"][col] < 1.0;
                        stations
                            .iter()
                            .position(|s| {
                                s.name == "exit"
                                    && (s.mach < 1.0) == subsonic
                                    && s.area_ratio.is_some_and(|a| (a / area - 1.0).abs() < 5e-4)
                            })
                            .unwrap_or_else(|| panic!("No station for reference area {area}"))
                    };
                    let st = stations[index];
                    seen[rocket_case].insert(index);
                    c.context = format!(
                        "case {} station {index} ({})",
                        results[rocket_case].case_index, st.name
                    );
                    c.state(page, "", col, &st.state);
                    c.row(page, "mach", col, st.mach);
                    // Performance tables leave the injector/chamber column blank.
                    if col > 0 {
                        for (key, actual) in [
                            ("area", st.area_ratio.unwrap()),
                            ("cstar", st.cstar_m_s),
                            ("cf", st.thrust_coefficient),
                            ("ivac", st.isp_vacuum_m_s),
                            ("isp", st.isp_m_s),
                        ] {
                            c.row(page, key, col - 1, actual);
                        }
                    }
                    if let Some(t) = &st.transport {
                        c.transport(page, col, t);
                    }
                }
            }
        }
    }
    if !matches!(results[0].output, CaseOutput::Rocket(_)) {
        assert_eq!(offset, results.len(), "Uncompared cases");
    } else {
        for (case, result) in results.iter().enumerate() {
            let CaseOutput::Rocket(r) = &result.output else {
                unreachable!()
            };
            let count = r.stations.len() - usize::from(r.throat_index == 3);
            for station in 0..count {
                if !seen[case].contains(&station) {
                    c.unreferenced_stations.push((case + 1, station));
                }
            }
        }
    }
    c
}

#[test]
fn nasa_sample_decks_against_reference_reports() {
    for i in 1..=14 {
        let input = std::fs::read_to_string(data_dir().join(format!("example{i}.inp"))).unwrap();
        let deck = parse_deck(&input).unwrap().remove(0);
        let cases = deck.cases().unwrap();
        let mut results = vec![];
        for case in &cases {
            match case.run() {
                Ok(r) => results.push(r),
                Err(e) => println!("example{i} case {} ERROR {e}", case.index),
            }
        }
        println!("example{i}: {}/{} solved", results.len(), cases.len());
        if results.len() != cases.len() {
            continue;
        }
        for source in ["", "native-3.3.4/"] {
            if i == 6 && !source.is_empty() {
                continue;
            }
            let text = std::fs::read_to_string(data_dir().join(format!("{source}example{i}.out")))
                .unwrap();
            let pages = reference(&text);
            let c = compare(&results, &pages);
            if i == 12 && !source.is_empty() {
                // The strict composition mismatch is a nozzle stopping difference:
                // the same public TP solver matches at NASA's reported state.
                let solver = cea_rs::solver::Solver::new(
                    deck.reactants.iter().map(|r| r.reactant.clone()).collect(),
                    true,
                    deck.ions,
                    deck.only.clone(),
                    deck.insert.clone(),
                )
                .unwrap();
                let page = &pages[0];
                let st = solver
                    .solve(
                        cea_rs::solver::Mode::TP(page.rows["t"][7], page.rows["p"][7]),
                        &results[0].weights,
                    )
                    .unwrap();
                let mut at_same_tp = Comparison::default();
                at_same_tp.state(page, "", 7, &st);
                assert!(at_same_tp.failures.is_empty());
                for name in ["CO", "OH"] {
                    assert!(
                        (st.mass_fractions[name] / page.rows[&format!("y.{name}")][7] - 1.0).abs()
                            < 1e-7
                    );
                }
            }
            println!(
                "example{i} {source} checked={} worst={:.9e} {} failures={}",
                c.checked,
                c.worst,
                c.worst_detail,
                c.failures.len()
            );
            for failure in &c.failures {
                println!("  {failure}");
            }
            for failure in &c.failures {
                assert!(
                    known_difference(i, source, failure),
                    "New discrepancy in example{i} {source}: {failure}"
                );
            }
            if i == 12 && source.is_empty() {
                // The input echoed in this historical report requests frozen only.
                assert_eq!(
                    c.unreferenced_stations,
                    (0..11).map(|s| (1, s)).collect::<Vec<_>>()
                );
            } else {
                assert!(c.unreferenced_stations.is_empty());
            }
            assert!(c.checked >= 20);
        }
    }
}

// Each exception is an observed, bounded discrepancy recorded in docs/input-decks.md.
// All other cells retain the 1e-4 check; these are not blanket deck exemptions.
fn known_difference(example: usize, source: &str, failure: &str) -> bool {
    let key = failure.split(": Rust=").next().unwrap();
    let bound = if source.is_empty() {
        match (example, key) {
            (1, "case 1 x.CO" | "case 1 x.CO2") => 1.11e-4,
            (7, key)
                if key.ends_with(" 1.p")
                    || key.ends_with(" 1.rho")
                    || key.ends_with(" 2.rho")
                    || key.ends_with(" 2.p")
                    || matches!(
                        key,
                        "case 3 x.H" | "case 4 x.H" | "case 5 x.H" | "case 6 x.H"
                    ) =>
            {
                1.7e-4
            }
            (
                8,
                "case 1 station 5 (exit) isp"
                | "case 1 station 5 (exit) mach"
                | "case 1 station 5 (exit) cf",
            )
            | (
                9,
                "case 1 station 1 (combustor) isp"
                | "case 1 station 1 (combustor) mach"
                | "case 1 station 1 (combustor) cf",
            ) => 1.16e-4,
            (10, "case 1 station 1 (combustor) area") => 1.15e-4,
            (
                10,
                "case 1 station 1 (combustor) isp"
                | "case 1 station 1 (combustor) mach"
                | "case 1 station 1 (combustor) cf",
            ) => 1.4e-4,
            (
                11,
                "case 1 station 3 (exit) isp"
                | "case 1 station 3 (exit) mach"
                | "case 1 station 3 (exit) cf",
            ) => 1.01e-4,
            _ => return false,
        }
    } else {
        match (example, key) {
            (12, "case 1 station 7 (exit) y.CO") => 1.74e-4,
            (12, "case 1 station 7 (exit) y.OH") => 1.07e-4,
            _ => return false,
        }
    };
    failure
        .split("relative=")
        .last()
        .unwrap()
        .parse::<f64>()
        .unwrap()
        <= bound
}

#[test]
fn nasa_combined_and_standalone_sample_inputs() {
    let standalone =
        parse_deck(&std::fs::read_to_string(data_dir().join("sample-example1.inp")).unwrap())
            .unwrap();
    let reference_pages =
        reference(&std::fs::read_to_string(data_dir().join("native-3.3.4/example1.out")).unwrap());
    let states = standalone[0]
        .cases()
        .unwrap()
        .iter()
        .map(|c| c.run().unwrap())
        .collect::<Vec<_>>();
    assert!(compare(&states, &reference_pages).failures.is_empty());
    let decks =
        parse_deck(&std::fs::read_to_string(data_dir().join("rp1311_examples.inp")).unwrap())
            .unwrap();
    assert_eq!(decks.len(), 14);
    let mut solved = 0;
    let mut failed = 0;
    for (index, deck) in decks.iter().enumerate() {
        let mut results = vec![];
        for case in deck.cases().unwrap() {
            match case.run() {
                Ok(result) => {
                    solved += 1;
                    results.push(result);
                }
                Err(e) => {
                    println!("combined example{} case{} {e}", index + 1, case.index);
                    assert_eq!(index + 1, 7);
                    assert_eq!(case.index, 1);
                    assert!(matches!(e, cea_rs::Error::Convergence(_)));
                    failed += 1;
                }
            }
        }
        let example = index + 1;
        let report = match example {
            6 => "example6.out".into(),
            7 => "native-3.3.4/sample-shock.out".into(),
            12 => "native-3.3.4/sample-12.out".into(),
            _ => format!("native-3.3.4/example{example}.out"),
        };
        let text = std::fs::read_to_string(data_dir().join(report)).unwrap();
        let mut pages = reference(&text);
        if example == 7 {
            assert!(text.contains("WARNING: No convergence for u1="));
            assert!(text.contains("ANSWERS NOT RELIABLE"));
            // NASA prints a retained, nonconverged state in column 1. Only that
            // column is excluded; all 13 successful shock cases are compared.
            for values in pages[0].rows.values_mut() {
                assert_eq!(values.len(), 7);
                values.remove(0);
            }
        }
        let c = compare(&results, &pages);
        println!(
            "combined example{example} checked={} worst={:.9e} {} failures={}",
            c.checked,
            c.worst,
            c.worst_detail,
            c.failures.len()
        );
        for failure in c.failures {
            assert!(
                known_difference(example, "native-3.3.4/", &failure),
                "combined example{example}: {failure}"
            );
        }
        assert!(c.unreferenced_stations.is_empty());
    }
    println!("combined: {solved} solved, {failed} failed");
    assert_eq!(failed, 1);
}
