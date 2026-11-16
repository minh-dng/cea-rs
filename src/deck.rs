//! NASA CEA input datasets, with SI values and no printed-report dependency.
//!
//! Derived from NASA CEA `source/{input,main,units}.f90` at
//! 2f79a647e85737742fb08a127da00d1ea1d069da (Apache-2.0); see NOTICE.
use crate::{Error, R, Result, thermo::Reactant};
use std::collections::BTreeMap;
use thermo_species::Species;
mod run;
pub use run::{Case, CaseInput, CaseOutput, DeckResult, run_deck};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Problem {
    TP,
    HP,
    SP,
    TV,
    UV,
    SV,
    Rocket,
    Shock,
    Detonation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Name,
    Fuel,
    Oxidant,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmountBasis {
    Mass,
    Moles,
}

#[derive(Debug, Clone)]
pub struct DeckReactant {
    pub reactant: Reactant,
    pub role: Role,
    pub amount: f64,
    pub basis: AmountBasis,
    pub temperature_k: Option<f64>,
    pub density_kg_m3: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ratio {
    OxidantFuel,
    FuelOxidant,
    PercentFuel,
    Phi,
    Chemical,
}

#[derive(Debug, Clone, Default)]
pub struct OutputOptions {
    pub transport: bool,
    pub trace: f64,
    pub mass_fractions: bool,
    /// Report preferences are metadata; returned numbers always use SI.
    pub calories: bool,
    pub long: bool,
    pub short: bool,
    pub debug: Vec<f64>,
    pub plot: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct RocketInput {
    pub finite_area: bool,
    pub ac_at: Option<f64>,
    pub mass_flux_kg_s_m2: Option<f64>,
    pub freeze_station: usize,
    pub temperature_estimate_k: Option<f64>,
    pub pressure_ratios: Vec<f64>,
    pub subsonic_area_ratios: Vec<f64>,
    pub supersonic_area_ratios: Vec<f64>,
}
impl Default for RocketInput {
    fn default() -> Self {
        Self {
            finite_area: false,
            ac_at: None,
            mass_flux_kg_s_m2: None,
            freeze_station: 1,
            temperature_estimate_k: None,
            pressure_ratios: vec![],
            subsonic_area_ratios: vec![],
            supersonic_area_ratios: vec![],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Deck {
    pub name: String,
    pub problem: Problem,
    pub reactants: Vec<DeckReactant>,
    pub pressures_pa: Vec<f64>,
    pub temperatures_k: Vec<f64>,
    pub volumes_m3_kg: Vec<f64>,
    pub enthalpies_j_kg: Vec<f64>,
    pub energies_j_kg: Vec<f64>,
    pub entropies_j_kg_k: Vec<f64>,
    pub ratios: Option<(Ratio, Vec<f64>)>,
    pub equilibrium: bool,
    pub frozen: bool,
    pub ions: bool,
    pub only: Option<Vec<String>>,
    pub omit: Vec<String>,
    pub insert: Vec<String>,
    pub rocket: RocketInput,
    pub shock_speeds: Vec<crate::shock::ShockSpeed>,
    pub incident: bool,
    pub reflected: bool,
    pub output: OutputOptions,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::Input(message.into())
}

fn unknown(dataset: &str, keyword: &str) -> Error {
    invalid(format!("{dataset}: unrecognized keyword {keyword:?}"))
}

fn number(word: &str) -> Option<f64> {
    word.replace(['d', 'D'], "e").parse().ok()
}

// Commas separate numeric lists, but belong to names such as C8H18(L),n-octa.
fn tokens(text: &str) -> Vec<String> {
    let mut words = vec![];
    let mut word = String::new();
    for c in text.chars().chain([' ']) {
        if c.is_whitespace() || c == '=' || (c == ',' && number(&word).is_some()) {
            if !word.is_empty() {
                words.push(word.trim_matches(',').to_string());
                word.clear();
            }
        } else if c != ',' || !word.is_empty() {
            word.push(c);
        }
    }
    words.retain(|s| !s.is_empty());
    words
}

struct Words {
    words: std::vec::IntoIter<String>,
}
impl Words {
    fn new(text: &str) -> Self {
        Self {
            words: tokens(text).into_iter(),
        }
    }
    fn next(&mut self) -> Option<String> {
        self.words.next()
    }
    fn value(&mut self, key: &str) -> Result<String> {
        self.next()
            .ok_or_else(|| invalid(format!("Missing value for {key}")))
    }
    fn numbers(&mut self, key: &str) -> Result<Vec<f64>> {
        let mut values = vec![];
        while let Some(n) = self.words.as_slice().first().and_then(|w| number(w)) {
            self.next();
            if !n.is_finite() {
                return Err(invalid(format!("{key} requires finite values")));
            }
            values.push(n);
        }
        if values.is_empty() {
            return Err(invalid(format!("Missing numeric value for {key}")));
        }
        Ok(values)
    }
    fn scalar(&mut self, key: &str) -> Result<f64> {
        let values = self.numbers(key)?;
        if values.len() != 1 {
            return Err(invalid(format!("{key} requires exactly one value")));
        }
        Ok(values[0])
    }
}

fn key_units(key: &str) -> (String, String) {
    let key = key.to_ascii_lowercase().replace('(', ",").replace(')', "");
    match key.split_once(',') {
        Some((key, units)) => (key.into(), units.into()),
        None => (key, String::new()),
    }
}

fn temperature(value: f64, units: &str) -> Result<f64> {
    let t = match units {
        "" | "k" => value,
        "c" => value + 273.15,
        "r" => value / 1.8,
        "f" => (value - 32.0) / 1.8 + 273.15,
        _ => return Err(unknown("temperature unit", units)),
    };
    positive(t, "temperature")
}

fn temperature_units<'a>(key: &str, units: &'a str) -> &'a str {
    if !units.is_empty() {
        return units;
    }
    ["f", "c", "r"]
        .into_iter()
        .find(|unit| key[1..].contains(unit))
        .unwrap_or("k")
}

fn positive(value: f64, key: &str) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(invalid(format!("{key} must be positive and finite")));
    }
    Ok(value)
}

fn pressure(value: f64, units: &str) -> Result<f64> {
    let factor = match units {
        "" | "bar" => 1e5,
        "atm" => 101325.0,
        "psi" | "psia" => 6894.75729,
        // NASA units.f90 deliberately uses this rounded conversion.
        "mmhg" => 133.3,
        _ => return Err(unknown("pressure unit", units)),
    };
    positive(value * factor, "pressure")
}

fn density(value: f64, units: &str) -> Result<f64> {
    let factor = match units {
        "" | "g/cc" | "g/cm**3" => 1000.0,
        "kg/m**3" => 1.0,
        _ => return Err(unknown("density unit", units)),
    };
    positive(value * factor, "density")
}

fn molar_energy(value: f64, units: &str) -> Result<f64> {
    let factor = match units {
        "" | "j" | "j/mol" | "j/mole" => 1.0,
        "kj" | "kj/mol" | "kj/mole" => 1000.0,
        "cal" | "cal/mol" | "cal/mole" => 4.184,
        "kcal" | "kcal/mol" | "kcal/mole" => 4184.0,
        _ => return Err(unknown("reactant energy unit", units)),
    };
    let value = value * factor;
    if !value.is_finite() {
        return Err(invalid("Nonfinite reactant enthalpy"));
    }
    Ok(value)
}

fn set_schedule(target: &mut Vec<f64>, values: Vec<f64>, key: &str) -> Result<()> {
    if !target.is_empty() {
        return Err(invalid(format!("Duplicate schedule {key}")));
    }
    *target = values;
    Ok(())
}

#[derive(Default)]
struct Builder {
    deck: Option<Deck>,
    reactants: Option<Vec<DeckReactant>>,
}

impl Builder {
    fn finish(self) -> Result<Deck> {
        let mut d = self
            .deck
            .ok_or_else(|| invalid("Missing problem dataset"))?;
        d.reactants = self
            .reactants
            .ok_or_else(|| invalid("Missing reactants dataset"))?;
        if d.only.is_some() && !d.omit.is_empty() {
            return Err(invalid("only and omit cannot be combined"));
        }
        if d.ratios.is_some() && d.reactants[0].role == Role::Name {
            return Err(invalid("Mixture ratio requires fuel and oxidant reactants"));
        }
        for (given, allowed, keyword) in [
            (
                !d.temperatures_k.is_empty(),
                matches!(
                    d.problem,
                    Problem::TP
                        | Problem::TV
                        | Problem::Rocket
                        | Problem::Shock
                        | Problem::Detonation
                ),
                "t",
            ),
            (
                !d.enthalpies_j_kg.is_empty(),
                matches!(d.problem, Problem::HP | Problem::Rocket),
                "h/r",
            ),
            (!d.energies_j_kg.is_empty(), d.problem == Problem::UV, "u/r"),
            (
                !d.entropies_j_kg_k.is_empty(),
                matches!(d.problem, Problem::SP | Problem::SV),
                "s/r",
            ),
            (
                !d.volumes_m3_kg.is_empty(),
                matches!(d.problem, Problem::TV | Problem::UV | Problem::SV),
                "v/rho",
            ),
            (
                !d.shock_speeds.is_empty() || d.incident || d.reflected,
                d.problem == Problem::Shock,
                "u1/mach1/incd/refl",
            ),
            (
                d.rocket.finite_area
                    || d.rocket.ac_at.is_some()
                    || d.rocket.mass_flux_kg_s_m2.is_some()
                    || d.rocket.temperature_estimate_k.is_some()
                    || !d.rocket.pressure_ratios.is_empty()
                    || !d.rocket.subsonic_area_ratios.is_empty()
                    || !d.rocket.supersonic_area_ratios.is_empty()
                    || d.rocket.freeze_station != 1,
                d.problem == Problem::Rocket,
                "fac/acat/ma/tcest/pi/p/subar/supar/nfz",
            ),
        ] {
            if given && !allowed {
                return Err(invalid(format!(
                    "Keyword {keyword} is not used by {:?}",
                    d.problem
                )));
            }
        }
        let volume = matches!(d.problem, Problem::TV | Problem::UV | Problem::SV);
        if volume && !d.pressures_pa.is_empty() {
            return Err(invalid("Keyword p is not used by a volume problem"));
        }
        if (volume && d.volumes_m3_kg.is_empty()) || (!volume && d.pressures_pa.is_empty()) {
            return Err(invalid(if volume {
                "Missing v or rho schedule"
            } else {
                "Missing p schedule"
            }));
        }
        if matches!(d.problem, Problem::TP | Problem::TV | Problem::Detonation)
            && d.temperatures_k.is_empty()
        {
            return Err(invalid("Missing t schedule"));
        }
        if matches!(d.problem, Problem::SP | Problem::SV) && d.entropies_j_kg_k.is_empty() {
            return Err(invalid("Missing s/r schedule"));
        }
        if d.problem == Problem::UV && d.energies_j_kg.is_empty() {
            return Err(invalid("Missing u/r schedule"));
        }
        if d.problem == Problem::Shock {
            if d.shock_speeds.is_empty() {
                return Err(invalid("Missing u1 or mach1 schedule"));
            }
            if d.pressures_pa.len() > d.shock_speeds.len()
                || d.temperatures_k.len() > d.shock_speeds.len()
            {
                return Err(invalid(
                    "Shock p/t schedules must not be longer than u1/mach1",
                ));
            }
        }
        if !(0.0..1.0).contains(&d.output.trace) {
            return Err(invalid("trace must be in [0, 1)"));
        }
        Ok(d)
    }
}

fn parse_problem(text: &str) -> Result<Deck> {
    let mut d = Deck {
        name: String::new(),
        problem: Problem::TP,
        reactants: vec![],
        pressures_pa: vec![],
        temperatures_k: vec![],
        volumes_m3_kg: vec![],
        enthalpies_j_kg: vec![],
        energies_j_kg: vec![],
        entropies_j_kg_k: vec![],
        ratios: None,
        equilibrium: false,
        frozen: false,
        ions: false,
        only: None,
        omit: vec![],
        insert: vec![],
        rocket: RocketInput::default(),
        shock_speeds: vec![],
        incident: false,
        reflected: false,
        output: OutputOptions::default(),
    };
    let mut problem = None;
    let mut words = Words::new(text);
    while let Some(keyword) = words.next() {
        let (literal, units) = key_units(&keyword);
        // NASA input.f90 tries fixed prefixes in descending length.
        let key = [
            "ions", "tces", "case", "mach", "mdot", "det", "f/o", "f/a", "fac", "iac", "h/r",
            "inc", "nfz", "nfr", "o/f", "phi", "ref", "rho", "rkt", "sub", "sup", "s/r", "u/r",
            "%f", "ac", "eq", "fr", "fz", "hp", "ph", "pc", "pi", "ro", "sh", "sp", "ps", "sv",
            "vs", "tp", "pt", "tv", "vt", "uv", "vu", "u1", "ma", "t", "p", "v", "r",
        ]
        .into_iter()
        .find(|prefix| literal.starts_with(prefix))
        .ok_or_else(|| unknown("problem", &keyword))?;
        if !units.is_empty()
            && !matches!(
                key,
                "p" | "t"
                    | "rho"
                    | "v"
                    | "o/f"
                    | "f/o"
                    | "f/a"
                    | "%f"
                    | "phi"
                    | "r"
                    | "sub"
                    | "sup"
                    | "ma"
                    | "mdot"
                    | "u1"
                    | "mach"
            )
        {
            return Err(unknown("problem", &keyword));
        }
        if matches!(key, "sub" | "sup") && !matches!(units.as_str(), "" | "ae/at") {
            return Err(unknown("area ratio unit", &units));
        }
        let kind = match key {
            "tp" | "pt" => Some(Problem::TP),
            "hp" | "ph" => Some(Problem::HP),
            "sp" | "ps" => Some(Problem::SP),
            "tv" | "vt" => Some(Problem::TV),
            "uv" | "vu" => Some(Problem::UV),
            "sv" | "vs" => Some(Problem::SV),
            "rkt" | "ro" => Some(Problem::Rocket),
            "sh" => Some(Problem::Shock),
            "det" => Some(Problem::Detonation),
            _ => None,
        };
        if let Some(kind) = kind {
            if problem.replace(kind).is_some() {
                return Err(invalid("Multiple problem types"));
            }
            continue;
        }
        match key {
            "case" => d.name = words.value(key)?,
            "ions" => d.ions = true,
            "eq" => d.equilibrium = true,
            "fr" | "fz" => d.frozen = true,
            "iac" => d.rocket.finite_area = false,
            "fac" => d.rocket.finite_area = true,
            "inc" => d.incident = true,
            "ref" => d.reflected = true,
            "p" => {
                let unit = if !units.is_empty() {
                    units.as_str()
                } else if literal[1..].contains("mmh") {
                    "mmhg"
                } else if literal[1..].contains("psi") {
                    "psi"
                } else if literal[1..].contains("atm") {
                    "atm"
                } else {
                    "bar"
                };
                let v = words
                    .numbers(&keyword)?
                    .into_iter()
                    .map(|n| pressure(n, unit))
                    .collect::<Result<_>>()?;
                set_schedule(&mut d.pressures_pa, v, &keyword)?;
            }
            "t" => {
                let unit = temperature_units(&literal, &units);
                let v = words
                    .numbers(&keyword)?
                    .into_iter()
                    .map(|n| temperature(n, unit))
                    .collect::<Result<_>>()?;
                set_schedule(&mut d.temperatures_k, v, &keyword)?;
            }
            "rho" | "v" => {
                let v = words
                    .numbers(&keyword)?
                    .into_iter()
                    .map(|n| {
                        if key == "rho" {
                            return Ok(1.0 / density(n, &units)?);
                        }
                        let factor = match units.as_str() {
                            "cc/g" | "cm**3/g" => 1e-3,
                            "m**3/kg" => 1.0,
                            _ => return Err(unknown("volume unit", &units)),
                        };
                        positive(n * factor, "volume")
                    })
                    .collect::<Result<_>>()?;
                set_schedule(&mut d.volumes_m3_kg, v, &keyword)?;
            }
            "h/r" | "u/r" | "s/r" => {
                if !units.is_empty() {
                    return Err(unknown("problem", &keyword));
                }
                let v = words
                    .numbers(&keyword)?
                    .into_iter()
                    .map(|n| n * R)
                    .collect::<Vec<_>>();
                if v.iter().any(|n| !n.is_finite()) {
                    return Err(invalid(format!("Nonfinite {key}")));
                }
                let target = match key {
                    "h/r" => &mut d.enthalpies_j_kg,
                    "u/r" => &mut d.energies_j_kg,
                    _ => &mut d.entropies_j_kg_k,
                };
                set_schedule(target, v, &keyword)?;
            }
            "o/f" | "f/o" | "f/a" | "%f" | "phi" | "r" => {
                if !matches!(units.as_str(), "" | "e" | "eq" | "eq.ratio") {
                    return Err(unknown("mixture ratio unit", &units));
                }
                let ratio = match key {
                    "o/f" => Ratio::OxidantFuel,
                    "f/o" | "f/a" => Ratio::FuelOxidant,
                    "%f" => Ratio::PercentFuel,
                    "phi" => Ratio::Phi,
                    _ => Ratio::Chemical,
                };
                if d.ratios
                    .replace((ratio, words.numbers(&keyword)?))
                    .is_some()
                {
                    return Err(invalid("Multiple mixture ratio schedules"));
                }
            }
            "pi" | "pc" => {
                set_schedule(
                    &mut d.rocket.pressure_ratios,
                    words.numbers(&keyword)?,
                    &keyword,
                )?;
            }
            "sub" => {
                set_schedule(
                    &mut d.rocket.subsonic_area_ratios,
                    words.numbers(&keyword)?,
                    &keyword,
                )?;
            }
            "sup" => {
                set_schedule(
                    &mut d.rocket.supersonic_area_ratios,
                    words.numbers(&keyword)?,
                    &keyword,
                )?;
            }
            "ac" => d.rocket.ac_at = Some(words.scalar(&keyword)?),
            "ma" | "mdot" => {
                if !matches!(units.as_str(), "" | "kg/s/m**2") {
                    return Err(unknown("mass flux unit", &units));
                }
                d.rocket.mass_flux_kg_s_m2 = Some(words.scalar(&keyword)?);
            }
            "nfz" | "nfr" => {
                let n = positive(words.scalar(&keyword)?, &keyword)?;
                if n.fract() != 0.0 || n > usize::MAX as f64 {
                    return Err(invalid("nfz must be a positive integer"));
                }
                d.rocket.freeze_station = n as usize;
            }
            "tces" => {
                // Pinned parse_prob reads this scalar directly, bypassing unit conversion.
                d.rocket.temperature_estimate_k =
                    Some(positive(words.scalar(&keyword)?, "temperature")?)
            }
            "u1" | "mach" => {
                if !d.shock_speeds.is_empty() {
                    return Err(invalid("Specify only one u1/mach1 schedule"));
                }
                if !units.is_empty() && !(key == "u1" && units == "m/s") {
                    return Err(unknown("shock speed unit", &units));
                }
                d.shock_speeds = words
                    .numbers(&keyword)?
                    .into_iter()
                    .map(|n| {
                        if key == "u1" {
                            crate::shock::ShockSpeed::Velocity(n)
                        } else {
                            crate::shock::ShockSpeed::Mach(n)
                        }
                    })
                    .collect();
            }
            _ => return Err(unknown("problem", &keyword)),
        }
    }
    d.problem = problem.ok_or_else(|| invalid("Missing problem type"))?;
    Ok(d)
}

#[derive(Default)]
struct ReactantInput {
    name: String,
    role: Option<Role>,
    amount: Option<(f64, AmountBasis)>,
    temperature: Option<f64>,
    enthalpy: Option<f64>,
    density: Option<f64>,
    molecular_weight: Option<f64>,
    formula: BTreeMap<String, f64>,
}

fn formula_mass(formula: &BTreeMap<String, f64>) -> Result<f64> {
    formula
        .iter()
        .map(|(symbol, n)| {
            // Neutral monatomic records carry the same atomic weights as NASA's table.
            let atomic = Species::all()
                .iter()
                .find(|sp| {
                    let f = sp.data().constituents();
                    f.len() == 1 && f[0].0 == 1.0 && format!("{:?}", f[0].1) == *symbol
                })
                .ok_or_else(|| invalid(format!("Unknown formula element {symbol}")))?;
            Ok(n * atomic.data().mw())
        })
        .sum()
}

impl ReactantInput {
    fn finish(self, basis: AmountBasis, needs_enthalpy: bool) -> Result<DeckReactant> {
        let named = Reactant::named(&self.name);
        if named.is_err() && needs_enthalpy && self.enthalpy.is_none() {
            return Err(invalid(format!(
                "Custom reactant {} requires h for this problem",
                self.name
            )));
        }
        let reactant = if self.formula.is_empty()
            && self.enthalpy.is_none()
            && self.molecular_weight.is_none()
        {
            named?
        } else {
            let formula_override = !self.formula.is_empty();
            let formula = if !formula_override {
                named
                    .as_ref()
                    .map_err(Clone::clone)?
                    .formula
                    .iter()
                    .map(|(n, c)| (format!("{c:?}"), *n))
                    .collect()
            } else {
                self.formula
            };
            let mw = match self.molecular_weight {
                Some(mw) => mw,
                None if !formula_override => named.as_ref().map_err(Clone::clone)?.molecular_weight,
                None => formula_mass(&formula)?,
            };
            let t = self
                .temperature
                .unwrap_or_else(|| named.as_ref().map_or(298.15, |r| r.reference_temperature));
            let h = match self.enthalpy {
                Some(h) => h,
                None if !needs_enthalpy => 0.0,
                None => named
                    .as_ref()
                    .ok()
                    .map(|r| r.enthalpy(t))
                    .transpose()?
                    .unwrap_or(0.0),
            };
            Reactant::custom(self.name, formula, mw, h, t)?
        };
        Ok(DeckReactant {
            reactant,
            role: self.role.expect("reactant role set by parser"),
            amount: self.amount.map_or(1.0, |a| a.0),
            basis,
            temperature_k: self.temperature,
            density_kg_m3: self.density,
        })
    }
}

fn parse_reactants(text: &str, needs_enthalpy: bool) -> Result<Vec<DeckReactant>> {
    let mut words = Words::new(text);
    let mut inputs = vec![];
    let mut current: Option<ReactantInput> = None;
    let mut basis = None;
    while let Some(keyword) = words.next() {
        let (key, units) = key_units(&keyword);
        let role = match key.get(..2) {
            Some("fu") => Some(Role::Fuel),
            Some("ox") => Some(Role::Oxidant),
            Some("na") => Some(Role::Name),
            _ => None,
        };
        if let Some(role) = role {
            if let Some(input) = current.take() {
                inputs.push(input);
            }
            current = Some(ReactantInput {
                name: words.value(&keyword)?,
                role: Some(role),
                ..ReactantInput::default()
            });
            continue;
        }
        let r = current
            .as_mut()
            .ok_or_else(|| unknown("reactants (expected name/fuel/oxid)", &keyword))?;
        if keyword.len() <= 2
            && keyword.starts_with(|c: char| c.is_ascii_uppercase())
            && keyword.chars().all(|c| c.is_ascii_alphabetic())
        {
            let n = if words
                .words
                .as_slice()
                .first()
                .and_then(|w| number(w))
                .is_some()
            {
                words.scalar(&keyword)?
            } else {
                1.0
            };
            *r.formula.entry(keyword.to_ascii_uppercase()).or_default() += n;
            continue;
        }
        let field = if matches!(key.as_str(), "mw" | "mwt" | "molwt" | "wt/mol" | "wt/mole") {
            "mw"
        } else if key.starts_with("den") {
            "rho"
        } else {
            match key.get(..1) {
                Some("m" | "w") => "amount",
                Some("t") => "t",
                Some("h" | "u") => "h",
                Some("r") => "rho",
                _ => return Err(unknown("reactants", &keyword)),
            }
        };
        match field {
            "amount" => {
                if !matches!(units.as_str(), "" | "%") {
                    return Err(unknown("reactant amount unit", &units));
                }
                let b = if key.starts_with('m') {
                    AmountBasis::Moles
                } else {
                    AmountBasis::Mass
                };
                if basis.replace(b).is_some_and(|previous| previous != b) {
                    return Err(invalid("Cannot mix mass and mole reactant amounts"));
                }
                let n = words.scalar(&keyword)?;
                if n < 0.0 {
                    return Err(invalid("Reactant amount must be nonnegative"));
                }
                r.amount = Some((n, b));
            }
            "t" => {
                r.temperature = Some(temperature(
                    words.scalar(&keyword)?,
                    temperature_units(&key, &units),
                )?);
            }
            "h" => r.enthalpy = Some(molar_energy(words.scalar(&keyword)?, &units)?),
            "rho" => r.density = Some(density(words.scalar(&keyword)?, &units)?),
            "mw" => {
                let factor = match units.as_str() {
                    "" | "g/mol" | "g/mole" | "kg/kmol" => 1.0,
                    "kg/mol" | "kg/mole" => 1000.0,
                    _ => return Err(unknown("molecular weight unit", &units)),
                };
                r.molecular_weight = Some(positive(
                    words.scalar(&keyword)? * factor,
                    "molecular weight",
                )?);
            }
            _ => unreachable!(),
        }
    }
    if let Some(input) = current {
        inputs.push(input);
    }
    if inputs.is_empty() {
        return Err(invalid("Empty reactants dataset"));
    }
    let named = inputs[0].role == Some(Role::Name);
    if inputs.iter().any(|r| (r.role == Some(Role::Name)) != named) {
        return Err(invalid("Cannot mix name with fuel/oxid reactants"));
    }
    let out = inputs
        .into_iter()
        .map(|r| r.finish(basis.unwrap_or(AmountBasis::Mass), needs_enthalpy))
        .collect::<Result<Vec<_>>>()?;
    positive(out.iter().map(|r| r.amount).sum(), "Reactant amount sum")?;
    Ok(out)
}

fn parse_output(text: &str, o: &mut OutputOptions) -> Result<()> {
    let mut words = Words::new(text);
    while let Some(key) = words.next() {
        let literal = key.to_ascii_lowercase();
        let prefix = ["massf", "long", "tran", "trac", "deb", "dbg", "cal"]
            .into_iter()
            .find(|prefix| literal.starts_with(prefix))
            .unwrap_or(&literal);
        match prefix {
            "tran" => o.transport = true,
            "trac" => o.trace = words.scalar(&key)?,
            "massf" => o.mass_fractions = true,
            "siunits" | "siunit" => o.calories = false,
            "cal" => o.calories = true,
            "long" => o.long = true,
            "short" => o.short = true,
            "deb" | "dbg" => o.debug = words.numbers(&key)?,
            "plot" => {
                while let Some(column) = words.words.as_slice().first() {
                    if !matches!(
                        column.to_ascii_lowercase().as_str(),
                        "aeat" | "t" | "p" | "ivac" | "isp" | "mach" | "cf" | "cstar"
                    ) {
                        break;
                    }
                    o.plot.push(words.next().unwrap());
                }
                if o.plot.is_empty() {
                    return Err(invalid("plot requires a supported column"));
                }
            }
            _ => return Err(unknown("output", &key)),
        }
    }
    Ok(())
}

fn dataset(word: &str) -> Option<&'static str> {
    let word = word.to_ascii_lowercase();
    if word.starts_with("end") {
        return Some("end");
    }
    match word.get(..4)? {
        "prob" => Some("problem"),
        "reac" => Some("reactants"),
        "only" => Some("only"),
        "omit" => Some("omit"),
        "inse" => Some("insert"),
        "outp" => Some("output"),
        _ => None,
    }
}

fn apply_dataset(builder: &mut Builder, kind: &str, text: &str) -> Result<()> {
    match kind {
        "problem" => {
            if builder.deck.is_some() {
                return Err(invalid("Duplicate problem dataset"));
            }
            builder.deck = Some(parse_problem(text)?);
        }
        "reactants" => {
            if builder.reactants.is_some() {
                return Err(invalid("Duplicate reactants dataset"));
            }
            let d = builder
                .deck
                .as_ref()
                .ok_or_else(|| invalid("Missing problem dataset"))?;
            let needs_enthalpy = matches!(d.problem, Problem::HP | Problem::Rocket)
                && d.enthalpies_j_kg.is_empty()
                && (d.problem != Problem::Rocket || d.temperatures_k.is_empty());
            builder.reactants = Some(parse_reactants(text, needs_enthalpy)?);
        }
        _ => return Err(unknown("dataset", kind)),
    }
    Ok(())
}

/// Read one or more `end`-terminated decks. Dataset order is arbitrary;
/// species names retain their case and punctuation. Unknown keywords fail.
pub fn parse_deck(text: &str) -> Result<Vec<Deck>> {
    let mut decks = vec![];
    let mut datasets: Vec<(&str, String)> = vec![];
    for (line_index, raw) in text.lines().enumerate() {
        let line = raw.split(['#', '!']).next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let word = line.split_whitespace().next().unwrap();
        if let Some(kind) = dataset(word) {
            if kind == "end" {
                let mut builder = Builder::default();
                // Problem/reactants may follow output and species datasets.
                for required in ["problem", "reactants"] {
                    for (kind, data) in &datasets {
                        if *kind == required {
                            apply_dataset(&mut builder, kind, data)?;
                        }
                    }
                }
                let d = builder
                    .deck
                    .as_mut()
                    .ok_or_else(|| invalid("Missing problem dataset"))?;
                for (kind, data) in &datasets {
                    match *kind {
                        "output" => parse_output(data, &mut d.output)?,
                        "only" => {
                            if d.only.is_some() {
                                return Err(invalid("Duplicate only dataset"));
                            }
                            d.only = Some(tokens(data));
                        }
                        "omit" => d.omit.extend(tokens(data)),
                        "insert" => d.insert.extend(tokens(data)),
                        _ => (),
                    }
                }
                decks.push(builder.finish()?);
                datasets.clear();
            } else {
                datasets.push((kind, line[word.len()..].to_string()));
            }
        } else if let Some((_, data)) = datasets.last_mut() {
            data.push(' ');
            data.push_str(line);
        } else {
            return Err(invalid(format!(
                "line {}: unrecognized dataset keyword {word:?}",
                line_index + 1
            )));
        }
    }
    if !datasets.is_empty() {
        return Err(invalid("Missing end keyword"));
    }
    if decks.is_empty() {
        return Err(invalid("No input decks"));
    }
    Ok(decks)
}
