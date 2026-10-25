//! Validated interface around the attributed Gibbs/Newton implementation.
use crate::{
    Error, R, Result, equilibrium_core,
    thermo::{Feed, Reactant, check_species_temperature},
};
use std::collections::{BTreeMap, BTreeSet};
use std::str::FromStr;
use thermo_species::{Constituent, Species};

/// All state quantities are SI: K, Pa, J/kg, J/(kg K), m³/kg.
#[derive(Debug, Clone, Copy)]
pub enum Mode {
    TP(f64, f64),
    HP(f64, f64),
    SP(f64, f64),
    TV(f64, f64),
    UV(f64, f64),
    SV(f64, f64),
}

#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct State {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub density_kg_m3: f64,
    pub enthalpy_j_kg: f64,
    pub energy_j_kg: f64,
    pub entropy_j_kg_k: f64,
    /// Valid only when equilibrium is true; frozen states use cp_frozen instead.
    pub cp_equilibrium_j_kg_k: f64,
    pub cp_frozen_j_kg_k: f64,
    pub gamma_s: f64,
    pub molecular_weight_kg_kmol: f64,
    pub iterations: usize,
    pub equilibrium: bool,
    pub element_residual: f64,
    pub constraint_residual: f64,
    /// NASA reporting convention: denominator includes gas and condensed moles.
    pub mole_fractions: BTreeMap<String, f64>,
    pub trace_extrapolated_species: Vec<String>,
    pub mass_fractions: BTreeMap<String, f64>,
    pub(crate) amounts: Vec<(f64, Species)>,
}
impl State {
    pub fn sound_speed(&self) -> f64 {
        (self.gamma_s * self.pressure_pa / self.density_kg_m3).sqrt()
    }
}

#[derive(Clone)]
pub struct Solver {
    pub reactants: Vec<Reactant>,
    pub condensed: bool,
    pub ions: bool,
    pub(crate) only: Option<Vec<Species>>,
    pub(crate) insert: Vec<Species>,
}
impl Solver {
    /// Construct a solver with an optional explicit product pool and phase seeds.
    pub fn new(
        reactants: Vec<Reactant>,
        condensed: bool,
        ions: bool,
        only: Option<Vec<String>>,
        insert: Vec<String>,
    ) -> Result<Self> {
        if reactants.is_empty() {
            return Err(Error::Input("At least one reactant is required".into()));
        }
        let parse = |name: &String| {
            Species::from_str(name).map_err(|_| Error::Input(format!("Unknown product {name}")))
        };
        let only = only
            .map(|v| v.iter().map(parse).collect::<Result<Vec<_>>>())
            .transpose()?;
        let insert = insert.iter().map(parse).collect::<Result<Vec<_>>>()?;
        if only.as_ref().is_some_and(Vec::is_empty) {
            return Err(Error::Input("Product pool must not be empty".into()));
        }
        for &sp in only.iter().flatten().chain(&insert) {
            if (!condensed && sp.data().phase() != 0)
                || (!ions
                    && sp
                        .data()
                        .constituents()
                        .iter()
                        .any(|(_, e)| *e == Constituent::E))
            {
                return Err(Error::Input(
                    "Product/seed conflicts with condensed or ion settings".into(),
                ));
            }
        }
        for &sp in &insert {
            if sp.data().phase() == 0 || only.as_ref().is_some_and(|v| !v.contains(&sp)) {
                return Err(Error::Input(
                    "Insert must contain condensed species in the product pool".into(),
                ));
            }
        }
        Ok(Self {
            reactants,
            condensed,
            ions,
            only,
            insert,
        })
    }

    pub(crate) fn feeds(&self, weights: &[f64]) -> Result<Vec<Feed>> {
        if weights.len() != self.reactants.len()
            || weights.iter().any(|w| !w.is_finite() || *w < 0.0)
        {
            return Err(Error::Input(
                "Weights must match reactants and be finite, nonnegative".into(),
            ));
        }
        let sum: f64 = weights.iter().sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(Error::Input(
                "Weights must have a finite positive sum".into(),
            ));
        }
        let feeds: Vec<_> = self
            .reactants
            .iter()
            .zip(weights)
            .filter(|(_, w)| **w > 0.0)
            .map(|(r, w)| Feed {
                amount: w / sum / r.molecular_weight,
                molecular_weight: r.molecular_weight,
                formula: r.formula.clone(),
            })
            .collect();
        let charge: f64 = feeds
            .iter()
            .map(|f| {
                f.amount
                    * f.formula
                        .iter()
                        .filter(|(_, c)| *c == Constituent::E)
                        .map(|(n, _)| n)
                        .sum::<f64>()
            })
            .sum();
        if charge.abs() > 1e-12
            || (!self.ions
                && feeds
                    .iter()
                    .any(|f| f.formula.iter().any(|(_, c)| *c == Constituent::E)))
        {
            return Err(Error::Input(
                "Reactants must be electrically neutral; ions require ions=True".into(),
            ));
        }
        let elements: BTreeSet<_> = feeds
            .iter()
            .flat_map(|f| f.formula.iter().map(|(_, c)| *c))
            .collect();
        let pool: Vec<_> =
            self.only
                .clone()
                .unwrap_or_else(|| Species::all().to_vec())
                .into_iter()
                .filter(|sp| {
                    (self.condensed || sp.data().phase() == 0)
                        && sp.data().constituents().iter().all(|(_, c)| {
                            elements.contains(c) || (*c == Constituent::E && self.ions)
                        })
                })
                .collect();
        if !pool.iter().any(|s| s.data().phase() == 0) {
            return Err(Error::Input(
                "Product pool has no compatible gas species".into(),
            ));
        }
        if let Some(requested) = &self.only {
            if requested.iter().any(|sp| !pool.contains(sp)) {
                return Err(Error::Input(
                    "Product formula contains an absent reactant element".into(),
                ));
            }
        }
        if self.insert.iter().any(|sp| !pool.contains(sp)) {
            return Err(Error::Input(
                "Insert contains an incompatible reactant element".into(),
            ));
        }
        for e in &elements {
            if !pool
                .iter()
                .any(|s| s.data().constituents().iter().any(|(_, c)| c == e))
            {
                return Err(Error::Input(format!("Product pool cannot conserve {e:?}")));
            }
        }
        Ok(feeds)
    }

    /// Compute reactant-stream enthalpy in J/kg with no temperature extrapolation.
    pub fn reactant_enthalpy(&self, weights: &[f64], temperatures: &[f64]) -> Result<f64> {
        self.feeds(weights)?;
        if temperatures.len() != weights.len() {
            return Err(Error::Input("Temperatures must match reactants".into()));
        }
        let total: f64 = weights.iter().sum();
        let mut h = 0.0;
        for ((r, w), t) in self.reactants.iter().zip(weights).zip(temperatures) {
            if *w > 0.0 {
                h += w / total / r.molecular_weight * r.properties(*t)?.1 * 1000.0;
            }
        }
        Ok(h)
    }

    /// Solve an equilibrium state, returning no result on failed verification.
    pub fn solve(&self, mode: Mode, weights: &[f64]) -> Result<State> {
        let feeds = self.feeds(weights)?;
        let (const_t, const_s, const_p, state1, state2) = match mode {
            Mode::TP(t, p) => (true, false, true, t, p / 1e5),
            Mode::HP(h, p) => (false, false, true, h / R, p / 1e5),
            Mode::SP(s, p) => (false, true, true, s / R, p / 1e5),
            Mode::TV(t, v) => (true, false, false, t, v),
            Mode::UV(u, v) => (false, false, false, u / R, v),
            Mode::SV(s, v) => (false, true, false, s / R, v),
        };
        if !state1.is_finite()
            || !state2.is_finite()
            || state2 <= 0.0
            || (const_t && !(300.0..=20000.0).contains(&state1))
        {
            return Err(Error::Domain(
                "Finite state inputs, positive pressure/volume and 300..20000 K are required"
                    .into(),
            ));
        }
        let d = equilibrium_core::solve_for_products_debug_seeded(
            &feeds,
            const_t,
            const_s,
            const_p,
            state1,
            state2,
            self.condensed,
            self.ions,
            self.only.as_deref(),
            Some(&self.insert),
        );
        if !d.converged {
            return Err(Error::Convergence(format!(
                "Gibbs/Newton solver failed after {} iterations",
                d.iterations_used
            )));
        }
        let mut amounts: Vec<_> =
            d.nj.iter()
                .copied()
                .zip(d.gas_species.iter().copied())
                .collect();
        amounts.extend(
            d.condensed_nj
                .iter()
                .copied()
                .zip(d.condensed_species.iter().copied())
                .filter(|(n, _)| *n > 0.0),
        );
        let mut out = State {
            temperature_k: d.temperature_k,
            pressure_pa: d.pressure_bar * 1e5,
            density_kg_m3: 1.0 / d.volume_m3_per_kg,
            enthalpy_j_kg: d.enthalpy_kj_per_kg * 1000.0,
            energy_j_kg: d.energy_kj_per_kg * 1000.0,
            entropy_j_kg_k: d.entropy_kj_per_kg_k * 1000.0,
            cp_equilibrium_j_kg_k: d.cp_eq_kj_per_kg_k * 1000.0,
            cp_frozen_j_kg_k: d.cp_frozen_kj_per_kg_k * 1000.0,
            gamma_s: d.gamma_s,
            molecular_weight_kg_kmol: d.molar_mass_kg_per_kmol,
            iterations: d.iterations_used,
            trace_extrapolated_species: vec![],
            equilibrium: true,
            element_residual: 0.0,
            constraint_residual: 0.0,
            mole_fractions: BTreeMap::new(),
            mass_fractions: BTreeMap::new(),
            amounts,
        };
        if !(300.0..=20000.0).contains(&out.temperature_k) {
            return Err(Error::Domain(
                "Solved temperature outside 300..20000 K".into(),
            ));
        }
        let positive = [
            out.temperature_k,
            out.pressure_pa,
            out.density_kg_m3,
            out.cp_equilibrium_j_kg_k,
            out.cp_frozen_j_kg_k,
            out.gamma_s,
            out.molecular_weight_kg_kmol,
        ];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || [out.enthalpy_j_kg, out.energy_j_kg, out.entropy_j_kg_k]
                .iter()
                .any(|v| !v.is_finite())
            || out.amounts.iter().any(|(n, _)| !n.is_finite() || *n < 0.0)
        {
            return Err(Error::Convergence(
                "Nonphysical or nonfinite converged state".into(),
            ));
        }
        // Reject any appreciable out-of-fit product instead of returning extrapolated properties.
        let n_total: f64 = out.amounts.iter().map(|(n, _)| n).sum();
        let mass: f64 = out.amounts.iter().map(|(n, s)| n * s.data().mw()).sum();
        for &(n, sp) in &out.amounts {
            if n / n_total > 1e-11 {
                check_species_temperature(sp, out.temperature_k)?;
            } else if n > 0.0 && check_species_temperature(sp, out.temperature_k).is_err() {
                out.trace_extrapolated_species
                    .push(sp.data().symbol().into());
            }
            out.mole_fractions
                .insert(sp.data().symbol().into(), n / n_total);
            out.mass_fractions
                .insert(sp.data().symbol().into(), n * sp.data().mw() / mass);
        }
        let elements: BTreeSet<_> = feeds
            .iter()
            .flat_map(|f| f.formula.iter().map(|(_, e)| *e))
            .collect();
        for e in elements.iter().chain(if self.ions {
            Some(&Constituent::E)
        } else {
            None
        }) {
            let initial: f64 = feeds
                .iter()
                .map(|f| {
                    f.amount
                        * f.formula
                            .iter()
                            .filter(|(_, c)| c == e)
                            .map(|(n, _)| n)
                            .sum::<f64>()
                })
                .sum();
            let final_n: f64 = out
                .amounts
                .iter()
                .map(|(n, sp)| {
                    n * sp
                        .data()
                        .constituents()
                        .iter()
                        .filter(|(_, c)| c == e)
                        .map(|(n, _)| n)
                        .sum::<f64>()
                })
                .sum();
            let residual = (initial - final_n).abs() / initial.abs().max(1e-8);
            out.element_residual = out.element_residual.max(residual);
        }
        let (actual, target) = match mode {
            Mode::TP(t, _) | Mode::TV(t, _) => (out.temperature_k, t),
            Mode::HP(h, _) => (out.enthalpy_j_kg, h),
            Mode::SP(s, _) | Mode::SV(s, _) => (out.entropy_j_kg_k, s),
            Mode::UV(u, _) => (out.energy_j_kg, u),
        };
        let scale = match mode {
            Mode::HP(..) | Mode::UV(..) => target
                .abs()
                .max(R * out.temperature_k / out.molecular_weight_kg_kmol)
                .max(1.0),
            _ => target.abs().max(1.0),
        };
        out.constraint_residual = (actual - target).abs() / scale;
        if out.element_residual > 2e-5 || out.constraint_residual > 2e-5 {
            return Err(Error::Convergence(format!(
                "Failed independent residual check: elements={}, constraint={}",
                out.element_residual, out.constraint_residual
            )));
        }
        Ok(out)
    }
}
