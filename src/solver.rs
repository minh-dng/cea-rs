//! Validated interface around the attributed Gibbs/Newton implementation.
//! NASA domain and reporting rules derive from Apache-2.0 NASA CEA; see provenance.
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
    /// NASA solve status; a nonconverged state retains NASA's available fields.
    pub converged: bool,
    /// NASA SUCCESS=0, NOT_CONVERGED=8, FORTRAN_ABORT=10.
    pub last_error: i32,
    pub(crate) abort_message: Option<String>,
    pub cv_equilibrium_j_kg_k: f64,
    pub cv_frozen_j_kg_k: f64,
    pub(crate) gas_moles: f64,
    pub(crate) post_processed: bool,
    pub(crate) transport_atomic: bool,
    pub(crate) dlnv_dlnt: f64,
    pub(crate) dlnv_dlnp: f64,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub density_kg_m3: f64,
    /// NASA can retain volume after clearing a failed shock's density.
    pub(crate) volume_m3_per_kg: f64,
    pub enthalpy_j_kg: f64,
    pub energy_j_kg: f64,
    pub gibbs_energy_j_kg: f64,
    pub entropy_j_kg_k: f64,
    /// Valid only when equilibrium is true; zero denotes NASA's phase-coexistence convention.
    pub cp_equilibrium_j_kg_k: f64,
    pub cp_frozen_j_kg_k: f64,
    pub gamma_s: f64,
    pub molecular_weight_kg_kmol: f64,
    pub(crate) mean_molecular_weight: f64,
    pub iterations: usize,
    pub equilibrium: bool,
    pub element_residual: f64,
    pub constraint_residual: f64,
    /// NASA reporting convention: denominator includes gas and condensed moles.
    pub mole_fractions: BTreeMap<String, f64>,
    pub trace_extrapolated_species: Vec<String>,
    /// All nonzero products evaluated outside their declared NASA fit intervals.
    pub extrapolated_species: Vec<String>,
    pub mass_fractions: BTreeMap<String, f64>,
    pub(crate) amounts: Vec<(f64, Species)>,
    pub(crate) equilibrium_seed: Option<equilibrium_core::EquilibriumSeed>,
}
impl State {
    pub(crate) fn iteration_enthalpy(&self) -> f64 {
        self.amounts.iter().fold(0.0, |sum, (n, s)| {
            n.mul_add(
                s.data().h_over_r(self.temperature_k) / self.temperature_k,
                sum,
            )
        }) * self.temperature_k
            * R
    }
    pub fn sound_speed(&self) -> f64 {
        (self.gamma_s * R * self.gas_moles * self.temperature_k).sqrt()
    }

    /// This composition frozen at `t` and `p`; its element balance carries over.
    pub(crate) fn frozen(&self, t: f64, p: f64) -> Result<State> {
        let mut st = frozen_state(&self.amounts, t, p, Some(self.gas_moles))?;
        st.molecular_weight_kg_kmol = self.molecular_weight_kg_kmol;
        st.mean_molecular_weight = self.mean_molecular_weight;
        st.mole_fractions.clone_from(&self.mole_fractions);
        st.mass_fractions.clone_from(&self.mass_fractions);
        st.element_residual = self.element_residual;
        Ok(st)
    }
}

/// Ideal-gas mixture with fixed species amounts (mol/kg) at `t` and `p`.
/// Frozen Cp is reported as both Cp values. `element_residual` is zero, which
/// suits amounts taken straight from a feed; use `State::frozen` otherwise.
/// Evaluates `t` and `p` directly, including nonphysical or nonfinite values.
pub(crate) fn frozen_state(
    amounts: &[(f64, Species)],
    t: f64,
    p: f64,
    gas_moles: Option<f64>,
) -> Result<State> {
    let total: f64 = amounts.iter().map(|(n, _)| n).sum();
    let ng: f64 = gas_moles.unwrap_or_else(|| {
        amounts
            .iter()
            .filter(|(_, sp)| sp.data().phase() == 0)
            .map(|(n, _)| n)
            .sum()
    });
    let mass: f64 = amounts.iter().map(|(n, sp)| n * sp.data().mw()).sum();
    let condensed: f64 = amounts
        .iter()
        .filter(|(_, sp)| sp.data().phase() != 0)
        .map(|(n, _)| n)
        .sum();
    let (mut h, mut s, mut cp) = (0.0, 0.0, 0.0);
    let (mut trace, mut extrapolated) = (vec![], vec![]);
    let mut x = BTreeMap::new();
    let mut y = BTreeMap::new();
    for &(n, sp) in amounts {
        // Absent species keep a zero entry, as in the state they were frozen from.
        x.insert(sp.data().symbol().into(), n / total);
        y.insert(sp.data().symbol().into(), n * sp.data().mw() / mass);
        if n <= 0.0 {
            continue;
        }
        // Same rule as Solver::solve: NASA extrapolates gas fits; report it.
        if check_species_temperature(sp, t).is_err() {
            extrapolated.push(sp.data().symbol().into());
            if n / total <= 1e-11 {
                trace.push(sp.data().symbol().into());
            }
        }
        let d = sp.data();
        h = n.mul_add(d.h_over_r(t) / t, h);
        cp = n.mul_add(d.cp_over_r(t), cp);
        let mixing = if d.phase() == 0 {
            (n / ng * p / 1e5).ln()
        } else {
            0.0
        };
        s = n.mul_add(d.s_over_r(t) - mixing, s);
    }
    h *= R * t;
    cp *= R;
    s *= R;
    let cv = cp - ng * R;
    Ok(State {
        converged: true,
        last_error: 0,
        abort_message: None,
        cv_equilibrium_j_kg_k: cv,
        cv_frozen_j_kg_k: cv,
        gas_moles: ng,
        post_processed: true,
        transport_atomic: false,
        dlnv_dlnt: 1.0,
        dlnv_dlnp: -1.0,
        temperature_k: t,
        pressure_pa: p,
        density_kg_m3: p / (ng * R * t),
        volume_m3_per_kg: ng * R * t / p,
        enthalpy_j_kg: h,
        energy_j_kg: h - ng * R * t,
        gibbs_energy_j_kg: h - t * s,
        entropy_j_kg_k: s,
        cp_equilibrium_j_kg_k: cp,
        cp_frozen_j_kg_k: cp,
        gamma_s: cp / cv,
        molecular_weight_kg_kmol: 1.0 / ng,
        mean_molecular_weight: (1.0 - condensed / total) / ng,
        iterations: 0,
        equilibrium: false,
        element_residual: 0.0,
        constraint_residual: 0.0,
        mole_fractions: x,
        trace_extrapolated_species: trace,
        extrapolated_species: extrapolated,
        mass_fractions: y,
        amounts: amounts.to_vec(),
        equilibrium_seed: None,
    })
}

#[derive(Clone)]
pub struct Solver {
    pub reactants: Vec<Reactant>,
    pub condensed: bool,
    pub ions: bool,
    pub(crate) only: Option<Vec<Species>>,
    pub(crate) insert: Vec<Species>,
    pub(crate) omit: Vec<Species>,
    pub(crate) trace: f64,
    pub(crate) transport: bool,
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
            omit: vec![],
            trace: 0.0,
            transport: false,
        })
    }

    /// Exclude named products from equilibrium.
    pub fn with_omit(mut self, names: Vec<String>) -> Result<Self> {
        // NASA filters product names exactly; names absent from the database are no-ops.
        self.omit = names
            .iter()
            .filter_map(|name| Species::from_str(name).ok())
            .collect();
        if self.insert.iter().any(|sp| self.omit.contains(sp)) {
            return Err(Error::Input(
                "An inserted phase cannot also be omitted".into(),
            ));
        }
        Ok(self)
    }

    /// Set NASA's trace convergence level; zero uses the default thresholds.
    pub fn with_trace(mut self, trace: f64) -> Result<Self> {
        if !trace.is_finite() || !(0.0..1.0).contains(&trace) {
            return Err(Error::Input("Trace must be finite and in [0, 1)".into()));
        }
        self.trace = trace;
        Ok(self)
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
                    !self.omit.contains(sp)
                        && (self.condensed || sp.data().phase() == 0)
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
        if let Some(requested) = &self.only
            && requested
                .iter()
                .any(|sp| !self.omit.contains(sp) && !pool.contains(sp))
        {
            return Err(Error::Input(
                "Product formula contains an absent reactant element".into(),
            ));
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
        if pool.len() == 1 {
            let mut amount: Option<f64> = None;
            for e in elements.iter().filter(|e| **e != Constituent::E) {
                let inventory: f64 = feeds
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
                let coefficient: f64 = pool[0]
                    .data()
                    .constituents()
                    .iter()
                    .filter(|(_, c)| c == e)
                    .map(|(n, _)| n)
                    .sum();
                let candidate = inventory / coefficient;
                if amount.is_some_and(|n| (candidate - n).abs() > 1e-10 * n.abs()) {
                    return Err(Error::Unsupported("Single-product pools with incompatible element inventories (native CEA can crash)".into()));
                }
                amount = Some(candidate);
            }
        }
        Ok(feeds)
    }

    /// Compute reactant-stream enthalpy in J/kg using NASA's fit/reference rule.
    pub fn reactant_enthalpy(&self, weights: &[f64], temperatures: &[f64]) -> Result<f64> {
        self.feeds(weights)?;
        if temperatures.len() != weights.len() {
            return Err(Error::Input("Temperatures must match reactants".into()));
        }
        let total: f64 = weights.iter().sum();
        let mut h = 0.0;
        for ((r, w), t) in self.reactants.iter().zip(weights).zip(temperatures) {
            if *w > 0.0 {
                h += w / total / r.molecular_weight * r.enthalpy(*t)? * 1000.0;
            }
        }
        if !h.is_finite() {
            return Err(Error::Domain("Nonfinite reactant enthalpy".into()));
        }
        Ok(h)
    }

    /// Solve an equilibrium state, retaining NASA's values and convergence status.
    pub fn solve(&self, mode: Mode, weights: &[f64]) -> Result<State> {
        self.solve_with_start(mode, weights, None)
    }

    pub(crate) fn initial_seed(
        &self,
        weights: &[f64],
        t: f64,
        amounts: Option<&[(f64, Species)]>,
    ) -> Result<equilibrium_core::EquilibriumSeed> {
        let feeds = self.feeds(weights)?;
        let elements: BTreeSet<_> = feeds
            .iter()
            .flat_map(|f| f.formula.iter().map(|(_, c)| *c))
            .collect();
        let mut species: Vec<_> =
            Species::all()
                .iter()
                .copied()
                .filter(|sp| {
                    (self.condensed || sp.data().phase() == 0)
                        && self.only.as_ref().is_none_or(|pool| pool.contains(sp))
                        && !self.omit.contains(sp)
                        && sp.data().constituents().iter().all(|(_, c)| {
                            elements.contains(c) || (self.ions && *c == Constituent::E)
                        })
                })
                .collect();
        species.sort_by_key(|s| (s.data().phase() != 0, s.data().symbol()));
        Ok(equilibrium_core::EquilibriumSeed::initial(
            species,
            t,
            amounts,
            &self.insert,
            self.trace,
        ))
    }

    pub(crate) fn solve_with_start(
        &self,
        mode: Mode,
        weights: &[f64],
        start: Option<&equilibrium_core::EquilibriumSeed>,
    ) -> Result<State> {
        let (first, second, fixed_temperature) = match mode {
            Mode::TP(t, p) | Mode::TV(t, p) => (t, p, true),
            Mode::HP(h, p) | Mode::UV(h, p) => (h, p, false),
            Mode::SP(s, p) | Mode::SV(s, p) => (s, p, false),
        };
        if !first.is_finite()
            || !second.is_finite()
            || second <= 0.0
            || (fixed_temperature && first <= 0.0)
        {
            return Err(Error::Domain(
                "Finite state inputs and positive temperature/pressure/volume are required".into(),
            ));
        }
        self.solve_state(mode, weights, start, false)
    }

    pub(crate) fn solve_state(
        &self,
        mode: Mode,
        weights: &[f64],
        start: Option<&equilibrium_core::EquilibriumSeed>,
        finalize_failure: bool,
    ) -> Result<State> {
        let feeds = self.feeds(weights)?;
        let (const_t, const_s, const_p, state1, state2) = match mode {
            Mode::TP(t, p) => (true, false, true, t, p / 1e5),
            Mode::HP(h, p) => (false, false, true, h / R, p / 1e5),
            Mode::SP(s, p) => (false, true, true, s / R, p / 1e5),
            Mode::TV(t, v) => (true, false, false, t, v),
            Mode::UV(u, v) => (false, false, false, u / R, v),
            Mode::SV(s, v) => (false, true, false, s / R, v),
        };
        let pool = self
            .only
            .clone()
            .unwrap_or_else(|| Species::all().to_vec())
            .into_iter()
            .filter(|sp| !self.omit.contains(sp))
            .collect::<Vec<_>>();
        let d = equilibrium_core::solve_for_products_debug_with_start(
            &feeds,
            const_t,
            const_s,
            const_p,
            state1,
            state2,
            self.condensed,
            self.ions,
            Some(&pool),
            Some(&self.insert),
            self.trace,
            start,
            finalize_failure,
        );
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
                .filter(|(n, _)| *n != 0.0),
        );
        let mut out = State {
            converged: d.converged,
            last_error: d.last_error,
            abort_message: d.abort_message.clone(),
            cv_equilibrium_j_kg_k: d.cv_eq_kj_per_kg_k * 1000.0,
            cv_frozen_j_kg_k: d.cv_frozen_kj_per_kg_k * 1000.0,
            gas_moles: d.n,
            post_processed: d.post_processed || finalize_failure,
            transport_atomic: false,
            dlnv_dlnt: d.dlnv_dlnt,
            dlnv_dlnp: d.dlnv_dlnp,
            temperature_k: d.temperature_k,
            pressure_pa: d.pressure_bar * 1e5,
            density_kg_m3: 1.0 / d.volume_m3_per_kg,
            volume_m3_per_kg: d.volume_m3_per_kg,
            enthalpy_j_kg: d.enthalpy_kj_per_kg * 1000.0,
            energy_j_kg: d.energy_kj_per_kg * 1000.0,
            gibbs_energy_j_kg: d.gibbs_kj_per_kg * 1000.0,
            entropy_j_kg_k: d.entropy_kj_per_kg_k * 1000.0,
            cp_equilibrium_j_kg_k: d.cp_eq_kj_per_kg_k * 1000.0,
            cp_frozen_j_kg_k: d.cp_frozen_kj_per_kg_k * 1000.0,
            gamma_s: d.gamma_s,
            molecular_weight_kg_kmol: d.molar_mass_kg_per_kmol,
            mean_molecular_weight: d.mean_molecular_weight_kg_per_kmol,
            iterations: d.iterations_used,
            trace_extrapolated_species: vec![],
            extrapolated_species: vec![],
            equilibrium: true,
            element_residual: 0.0,
            constraint_residual: 0.0,
            mole_fractions: BTreeMap::new(),
            mass_fractions: BTreeMap::new(),
            amounts,
            equilibrium_seed: Some(d.seed),
        };
        // NASA extrapolates gas fits and checks condensed-phase validity in the active set.
        let n_total: f64 = out.amounts.iter().map(|(n, _)| n).sum();
        let mass: f64 = out.amounts.iter().map(|(n, s)| n * s.data().mw()).sum();
        for &(n, sp) in &out.amounts {
            if n > 0.0 && check_species_temperature(sp, out.temperature_k).is_err() {
                out.extrapolated_species.push(sp.data().symbol().into());
                if n / n_total <= 1e-11 {
                    out.trace_extrapolated_species
                        .push(sp.data().symbol().into());
                }
            }
            out.mole_fractions
                .insert(sp.data().symbol().into(), n / n_total);
            out.mass_fractions
                .insert(sp.data().symbol().into(), n * sp.data().mw() / mass);
        }
        for sp in &d.condensed_species {
            out.mole_fractions
                .entry(sp.data().symbol().into())
                .or_insert(0.0 / n_total);
            out.mass_fractions
                .entry(sp.data().symbol().into())
                .or_insert(0.0 / mass);
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
        if !d.post_processed && !finalize_failure {
            let (t, gas_moles, amounts, seed, iterations) = (
                out.temperature_k,
                out.gas_moles,
                out.amounts,
                out.equilibrium_seed,
                out.iterations,
            );
            out = if let Some(previous) = start.and_then(|s| s.reported.as_deref()) {
                previous.clone()
            } else {
                State {
                    temperature_k: t,
                    pressure_pa: 0.0,
                    density_kg_m3: 0.0,
                    volume_m3_per_kg: 0.0,
                    enthalpy_j_kg: 0.0,
                    energy_j_kg: 0.0,
                    gibbs_energy_j_kg: 0.0,
                    entropy_j_kg_k: 0.0,
                    cp_equilibrium_j_kg_k: 0.0,
                    cp_frozen_j_kg_k: 0.0,
                    cv_equilibrium_j_kg_k: 0.0,
                    cv_frozen_j_kg_k: 0.0,
                    gamma_s: 0.0,
                    molecular_weight_kg_kmol: 0.0,
                    mean_molecular_weight: 0.0,
                    mole_fractions: out
                        .mole_fractions
                        .keys()
                        .map(|n| (n.clone(), 0.0))
                        .collect(),
                    mass_fractions: out
                        .mass_fractions
                        .keys()
                        .map(|n| (n.clone(), 0.0))
                        .collect(),
                    amounts: vec![],
                    equilibrium_seed: None,
                    ..out
                }
            };
            out.temperature_k = t;
            out.gas_moles = gas_moles;
            out.amounts = amounts;
            out.equilibrium_seed = seed;
            out.iterations = iterations;
            out.post_processed = false;
            out.converged = d.converged;
            out.last_error = d.last_error;
            out.abort_message = d.abort_message;
            out.equilibrium_seed.as_mut().unwrap().reported =
                start.and_then(|s| s.reported.clone());
        }
        if out.converged && out.post_processed {
            if self.transport {
                let tr = out.transport()?;
                out.cp_frozen_j_kg_k = tr.cp_frozen_j_kg_k;
                out.cv_frozen_j_kg_k = tr.cp_frozen_j_kg_k - out.gas_moles * R;
                out.equilibrium_seed.as_mut().unwrap().has_transport_seed |=
                    tr.viscosity_pa_s > 0.0;
            }
            let mut reported = out.clone();
            reported.equilibrium_seed = None;
            out.equilibrium_seed.as_mut().unwrap().reported = Some(Box::new(reported));
        }
        Ok(out)
    }
}
