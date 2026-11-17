//! Mixture conversions derived from NASA CEA mixture.f90 (Apache-2.0).
//! State properties use NASA's ideal-gas/zero-volume condensed convention.
//! Ratio conversions retain f64 precision; native cea==3.3.4 Python wrappers
//! first round scalar ratios through float32 (for example r=0.9 differs by
//! 2.65e-8 relative). Exactly representable scalar goldens agree to 1e-12.
//! Frozen Cv is Cp-R/M at the current state. Native frozen rocket stations
//! can retain an earlier Cv or report zero, so those diagnostics are not a
//! usable frozen-Cv oracle. Equilibrium Cv is reconstructed without core edits.
use crate::{
    Error, R, Result,
    solver::{Solver, State},
};
use std::collections::BTreeMap;

fn positive(value: f64, name: &str, allow_zero: bool) -> Result<f64> {
    if !value.is_finite() || value < 0.0 || (!allow_zero && value == 0.0) {
        return Err(Error::Input(format!(
            "{name} must be finite and {}",
            if allow_zero {
                "nonnegative"
            } else {
                "positive"
            }
        )));
    }
    Ok(value)
}

impl Solver {
    pub(crate) fn check_mixture_vector(&self, values: &[f64]) -> Result<f64> {
        if values.len() != self.reactants.len() || values.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(Error::Input(
                "Amounts must match reactants and be finite, nonnegative".into(),
            ));
        }
        positive(values.iter().sum(), "Amount sum", true)
    }

    /// Convert masses to mole amounts without normalization (kg to kmol).
    pub fn weights_to_moles(&self, weights: &[f64]) -> Result<Vec<f64>> {
        self.check_mixture_vector(weights)?;
        let moles: Vec<_> = weights
            .iter()
            .zip(&self.reactants)
            .map(|(w, r)| w / r.molecular_weight)
            .collect();
        self.check_mixture_vector(&moles)?;
        Ok(moles)
    }

    /// Convert mole amounts to masses without normalization (kmol to kg).
    pub fn moles_to_weights(&self, moles: &[f64]) -> Result<Vec<f64>> {
        self.check_mixture_vector(moles)?;
        let weights: Vec<_> = moles
            .iter()
            .zip(&self.reactants)
            .map(|(n, r)| n * r.molecular_weight)
            .collect();
        self.check_mixture_vector(&weights)?;
        Ok(weights)
    }

    /// NASA's scaled weights; normalize only when supplying mass fractions.
    /// Both stream vectors use the solver's complete reactant order.
    pub fn of_ratio_to_weights(
        &self,
        oxidant_weights: &[f64],
        fuel_weights: &[f64],
        of_ratio: f64,
    ) -> Result<Vec<f64>> {
        positive(of_ratio, "O/F", true)?;
        let ox = positive(
            self.check_mixture_vector(oxidant_weights)?,
            "Oxidant sum",
            false,
        )?;
        let fu = positive(self.check_mixture_vector(fuel_weights)?, "Fuel sum", false)?;
        let mass: f64 = self
            .reactants
            .iter()
            .zip(oxidant_weights.iter().zip(fuel_weights))
            .map(|(r, (o, f))| (o + f) * r.molecular_weight)
            .sum();
        let weights: Vec<_> = oxidant_weights
            .iter()
            .zip(fuel_weights)
            .map(|(o, f)| {
                mass * (f / fu / (1.0 + of_ratio) + o / ox * (of_ratio / (1.0 + of_ratio)))
            })
            .collect();
        self.check_mixture_vector(&weights)?;
        Ok(weights)
    }

    /// Negative and positive NASA valences, keyed by element (E has zero).
    pub fn get_valence(&self) -> Result<BTreeMap<String, (f64, f64)>> {
        self.reactants
            .iter()
            .flat_map(|r| &r.formula)
            .map(|(_, e)| {
                let name = format!("{e:?}");
                let v = atom_valence(&name)?;
                Ok((name, (v.min(0.0), v.max(0.0))))
            })
            .collect()
    }

    fn stream_valence(&self, weights: &[f64]) -> Result<(f64, f64)> {
        let sum = positive(self.check_mixture_vector(weights)?, "Stream sum", false)?;
        let valences = self.get_valence()?;
        let mut minus = 0.0;
        let mut plus = 0.0;
        for (r, w) in self.reactants.iter().zip(weights) {
            for (count, e) in &r.formula {
                let (vm, vp) = valences[&format!("{e:?}")];
                minus += count * w / sum / r.molecular_weight * vm;
                plus += count * w / sum / r.molecular_weight * vp;
            }
        }
        Ok((minus, plus))
    }

    pub fn chem_eq_ratio_to_of_ratio(
        &self,
        oxidant_weights: &[f64],
        fuel_weights: &[f64],
        chem_eq_ratio: f64,
    ) -> Result<f64> {
        positive(chem_eq_ratio, "Chemical equivalence ratio", true)?;
        let (om, op) = self.stream_valence(oxidant_weights)?;
        let (fm, fp) = self.stream_valence(fuel_weights)?;
        positive(
            -(fp + chem_eq_ratio * fm) / (op + chem_eq_ratio * om),
            "Resulting O/F",
            true,
        )
    }

    pub fn of_ratio_to_chem_eq_ratio(
        &self,
        oxidant_weights: &[f64],
        fuel_weights: &[f64],
        of_ratio: f64,
    ) -> Result<f64> {
        positive(of_ratio, "O/F", true)?;
        let (om, op) = self.stream_valence(oxidant_weights)?;
        let (fm, fp) = self.stream_valence(fuel_weights)?;
        positive(
            -(fp + op * of_ratio) / (fm + om * of_ratio),
            "Resulting chemical equivalence ratio",
            true,
        )
    }

    pub fn weight_eq_ratio_to_of_ratio(
        &self,
        oxidant_weights: &[f64],
        fuel_weights: &[f64],
        weight_eq_ratio: f64,
    ) -> Result<f64> {
        positive(weight_eq_ratio, "Phi", false)?;
        let (om, op) = self.stream_valence(oxidant_weights)?;
        let (fm, fp) = self.stream_valence(fuel_weights)?;
        positive(
            -(fm + fp) / (weight_eq_ratio * (om + op)),
            "Resulting O/F",
            true,
        )
    }

    pub fn of_ratio_to_weight_eq_ratio(
        &self,
        oxidant_weights: &[f64],
        fuel_weights: &[f64],
        of_ratio: f64,
    ) -> Result<f64> {
        positive(of_ratio, "O/F", false)?;
        self.weight_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, of_ratio)
    }
}

// NASA atomic_data.f90 valences; zero-valence gases and E do not oxidize fuel.
fn atom_valence(name: &str) -> Result<f64> {
    Ok(match name {
        "E" | "HE" | "NE" | "AR" | "KR" | "XE" | "RN" | "N" => 0.0,
        "H" | "IH" | "D" | "LI" | "NA" | "K" | "RB" | "AG" | "CS" | "TL" | "FR" => 1.0,
        "BE" | "MG" | "MN" | "CA" | "CO" | "NI" | "CU" | "ZN" | "SR" | "PD" | "CD" | "BA"
        | "HG" | "PB" | "PO" | "RA" => 2.0,
        "B" | "AL" | "SC" | "CR" | "FE" | "GA" | "AS" | "Y" | "RU" | "RH" | "IN" | "SB" | "LA"
        | "CE" | "PR" | "ND" | "PM" | "SM" | "EU" | "GD" | "TB" | "DY" | "HO" | "ER" | "TM"
        | "YB" | "LU" | "AU" | "BI" | "AC" | "AM" | "CM" | "BK" | "CF" | "ES" => 3.0,
        "C" | "IC" | "SI" | "S" | "TI" | "GE" | "SE" | "ZR" | "SN" | "TE" | "HF" | "OS" | "IR"
        | "PT" | "TH" | "PU" => 4.0,
        "P" | "V" | "NB" | "TA" | "PA" | "NP" => 5.0,
        "MO" | "W" | "U" => 6.0,
        "TC" | "RE" => 7.0,
        "O" | "IO" => -2.0,
        "F" | "CL" | "BR" | "I" | "AT" => -1.0,
        _ => return Err(Error::Data(format!("No NASA valence for {name}"))),
    })
}

/// Gibbs energy in J/kg.
pub fn gibbs_energy(state: &State) -> f64 {
    state.enthalpy_j_kg - state.temperature_k * state.entropy_j_kg_k
}

/// Frozen Cv in J/(kg K), including condensed heat capacity.
pub fn cv_frozen(state: &State) -> f64 {
    state.cp_frozen_j_kg_k - R / state.molecular_weight_kg_kmol
}

/// Equilibrium Cv in J/(kg K), reconstructed from the untruncated TP Jacobian.
pub fn cv_equilibrium(state: &State) -> Result<f64> {
    crate::derivatives::cv_equilibrium(state)
}

pub fn specific_volume(state: &State) -> f64 {
    1.0 / state.density_kg_m3
}

/// NASA M = 1/n_gas; includes condensed mass in the one-kg mixture basis.
pub fn molecular_weight_m(state: &State) -> f64 {
    state.molecular_weight_kg_kmol
}

/// NASA MW = M times the gas fraction of total gas plus condensed moles.
pub fn molecular_weight_mw(state: &State) -> f64 {
    let total: f64 = state.amounts.iter().map(|(n, _)| n).sum();
    let condensed: f64 = state
        .amounts
        .iter()
        .filter(|(_, s)| s.data().phase() != 0)
        .map(|(n, _)| n)
        .sum();
    state.molecular_weight_kg_kmol * (1.0 - condensed / total)
}
