//! Reactant properties and element bookkeeping in SI units.
use crate::{Error, R, Result};
use std::collections::BTreeMap;
use std::str::FromStr;
use std::sync::OnceLock;
use thermo_species::{Constituent, Species};

#[derive(Clone, Debug)]
pub struct Feed {
    pub amount: f64,
    pub molecular_weight: f64,
    pub formula: Vec<(f64, Constituent)>,
}
impl Feed {
    pub fn from_species(amount: f64, sp: Species) -> Self {
        Self {
            amount,
            molecular_weight: sp.data().mw(),
            formula: sp.data().constituents().to_vec(),
        }
    }
}

pub fn element_amounts(feeds: &[Feed], elements: &[Constituent]) -> Vec<f64> {
    let mass: f64 = feeds.iter().map(|f| f.amount * f.molecular_weight).sum();
    elements
        .iter()
        .map(|e| {
            feeds
                .iter()
                .map(|f| {
                    f.amount
                        * f.formula
                            .iter()
                            .filter(|(_, c)| c == e)
                            .map(|(n, _)| n)
                            .sum::<f64>()
                        / mass
                })
                .sum()
        })
        .collect()
}

fn constituent(symbol: &str) -> Result<Constituent> {
    let symbol = symbol.to_ascii_uppercase();
    Species::all()
        .iter()
        .flat_map(|s| s.data().constituents())
        .map(|(_, c)| *c)
        .find(|c| format!("{c:?}") == symbol)
        .ok_or_else(|| Error::Input(format!("Unknown element {symbol}")))
}

#[derive(Clone, Debug)]
struct Fit {
    lo: f64,
    hi: f64,
    a: [f64; 9],
}

#[cfg_attr(
    feature = "python",
    pyo3::pyclass(name = "Reactant", frozen, skip_from_py_object)
)]
#[derive(Clone, Debug)]
pub struct Reactant {
    pub name: String,
    pub molecular_weight: f64,
    pub formula: Vec<(f64, Constituent)>,
    pub reference_temperature: f64,
    reference_enthalpy: f64,
    fits: Vec<Fit>,
    species: Option<Species>,
}

impl Reactant {
    /// Load a NASA species, including reference-only liquid fuels.
    pub fn named(name: &str) -> Result<Self> {
        if let Ok(sp) = Species::from_str(name) {
            return Ok(Self {
                name: name.into(),
                molecular_weight: sp.data().mw(),
                formula: sp.data().constituents().to_vec(),
                reference_temperature: 298.15,
                reference_enthalpy: sp.data().h_formation(),
                fits: vec![],
                species: Some(sp),
            });
        }
        static EXTRA: OnceLock<std::result::Result<Vec<Reactant>, String>> = OnceLock::new();
        let extra = EXTRA.get_or_init(|| parse_reactants().map_err(|e| e.to_string()));
        extra
            .as_ref()
            .map_err(|e| Error::Data(e.clone()))?
            .iter()
            .find(|r| r.name == name)
            .cloned()
            .ok_or_else(|| Error::Input(format!("Unknown NASA reactant {name}")))
    }

    /// Define a custom reference enthalpy (J/mol), formula and molar mass (g/mol).
    pub fn custom(
        name: String,
        formula: BTreeMap<String, f64>,
        molecular_weight: f64,
        enthalpy_j_mol: f64,
        reference_temperature_k: f64,
    ) -> Result<Self> {
        if name.trim().is_empty()
            || formula.is_empty()
            || !molecular_weight.is_finite()
            || molecular_weight <= 0.0
            || !enthalpy_j_mol.is_finite()
            || !reference_temperature_k.is_finite()
            || reference_temperature_k <= 0.0
        {
            return Err(Error::Input("Custom reactant requires name, formula, positive mass/temperature and finite enthalpy".into()));
        }
        let mut counts = BTreeMap::new();
        for (symbol, n) in formula {
            let c = constituent(&symbol)?;
            if !n.is_finite() || n == 0.0 || (n < 0.0 && c != Constituent::E) {
                return Err(Error::Input(
                    "Element counts must be finite and positive (signed E denotes charge)".into(),
                ));
            }
            *counts.entry(c).or_insert(0.0) += n;
        }
        Ok(Self {
            name,
            molecular_weight,
            formula: counts.into_iter().map(|(c, n)| (n, c)).collect(),
            reference_temperature: reference_temperature_k,
            reference_enthalpy: enthalpy_j_mol,
            fits: vec![],
            species: None,
        })
    }

    /// Return molar heat capacity, enthalpy and standard entropy in J units.
    pub fn properties(&self, t: f64) -> Result<(Option<f64>, f64, Option<f64>)> {
        if !t.is_finite() || t <= 0.0 {
            return Err(Error::Domain(
                "Temperature must be positive and finite".into(),
            ));
        }
        if let Some(sp) = self.species {
            check_species_temperature(sp, t)?;
            let d = sp.data();
            return Ok((
                Some(d.cp_over_r(t) * R / 1000.0),
                d.h_over_r(t) * R / 1000.0,
                Some(d.s_over_r(t) * R / 1000.0),
            ));
        }
        if self.fits.is_empty() {
            if (t - self.reference_temperature).abs() > 1e-7 {
                return Err(Error::Domain(format!(
                    "{} has enthalpy only at {} K; no heat capacity is available",
                    self.name, self.reference_temperature
                )));
            }
            return Ok((None, self.reference_enthalpy, None));
        }
        let f = self
            .fits
            .iter()
            .find(|f| t >= f.lo && t <= f.hi)
            .ok_or_else(|| {
                Error::Domain(format!(
                    "{}: {t} K is outside the fitted intervals",
                    self.name
                ))
            })?;
        let a = f.a;
        let cp = a[0] / t.powi(2)
            + a[1] / t
            + a[2]
            + a[3] * t
            + a[4] * t * t
            + a[5] * t.powi(3)
            + a[6] * t.powi(4);
        let h = -a[0] / t
            + a[1] * t.ln()
            + a[2] * t
            + a[3] * t * t / 2.0
            + a[4] * t.powi(3) / 3.0
            + a[5] * t.powi(4) / 4.0
            + a[6] * t.powi(5) / 5.0
            + a[7];
        let s = -a[0] / (2.0 * t * t) - a[1] / t
            + a[2] * t.ln()
            + a[3] * t
            + a[4] * t * t / 2.0
            + a[5] * t.powi(3) / 3.0
            + a[6] * t.powi(4) / 4.0
            + a[8];
        Ok((Some(cp * R / 1000.0), h * R / 1000.0, Some(s * R / 1000.0)))
    }
}

pub(crate) fn check_species_temperature(sp: Species, t: f64) -> Result<()> {
    if !t.is_finite()
        || !sp
            .data()
            .temperature_data()
            .iter()
            .any(|(lo, hi, _)| t >= *lo - 1e-9 && t <= *hi + 1e-9)
    {
        return Err(Error::Domain(format!(
            "{}: {t} K is outside the NASA fitted intervals",
            sp.data().symbol()
        )));
    }
    Ok(())
}

fn number(line: &str, lo: usize, hi: usize) -> Result<f64> {
    line.get(lo..hi)
        .ok_or_else(|| Error::Data("Truncated NASA database record".into()))?
        .trim()
        .replace('D', "e")
        .parse::<f64>()
        .map_err(|e| Error::Data(e.to_string()))
}

fn parse_reactants() -> Result<Vec<Reactant>> {
    let text = include_str!("../vendor/thermo-species/data/thermo.inp");
    let mut lines = text.lines().skip_while(|l| *l != "END PRODUCTS").skip(1);
    let mut out = vec![];
    while let Some(header) = lines.next() {
        if header == "END REACTANTS" {
            break;
        }
        let meta = lines
            .next()
            .ok_or_else(|| Error::Data("Missing NASA metadata".into()))?;
        let count = number(meta, 0, 2)? as usize;
        let mut formula = BTreeMap::new();
        for i in 0..5 {
            let start = 10 + 8 * i;
            let symbol = meta[start..start + 2].trim();
            let n = number(meta, start + 2, start + 8)?;
            if !symbol.is_empty() && n != 0.0 {
                formula.insert(symbol.to_string(), n);
            }
        }
        let mut r = Reactant::custom(
            header[..18].trim().into(),
            formula,
            number(meta, 52, 65)?,
            number(meta, 65, 80)?,
            298.15,
        )?;
        if count == 0 {
            r.reference_temperature = number(
                lines
                    .next()
                    .ok_or_else(|| Error::Data("Missing NASA reference state".into()))?,
                0,
                11,
            )?;
        } else {
            for _ in 0..count {
                let range = lines
                    .next()
                    .ok_or_else(|| Error::Data("Missing NASA range".into()))?;
                let row1 = lines
                    .next()
                    .ok_or_else(|| Error::Data("Missing NASA coefficients".into()))?;
                let row2 = lines
                    .next()
                    .ok_or_else(|| Error::Data("Missing NASA coefficients".into()))?;
                let mut a = [0.0; 9];
                for (i, x) in a[..5].iter_mut().enumerate() {
                    *x = number(row1, 16 * i, 16 * (i + 1))?;
                }
                // The CEA format reserves column 3; reject a different layout.
                let reserved = row2
                    .get(32..48)
                    .ok_or_else(|| Error::Data("Truncated NASA coefficient row".into()))?
                    .trim();
                if !reserved.is_empty() && number(row2, 32, 48)? != 0.0 {
                    return Err(Error::Data(
                        "Nonzero reserved NASA coefficient column".into(),
                    ));
                }
                for (i, j) in [0, 1, 3, 4].into_iter().enumerate() {
                    a[5 + i] = number(row2, 16 * j, 16 * (j + 1))?;
                }
                r.fits.push(Fit {
                    lo: number(range, 0, 11)?,
                    hi: number(range, 11, 22)?,
                    a,
                });
            }
        }
        out.push(r);
    }
    Ok(out)
}
