//! Equilibrium sensitivities derived from NASA CEA equilibrium.f90 (Apache-2.0).
//!
//! Every response vector is ordered [state1, state2, reactant weights...].
//! Inputs and responses use SI; mole amounts are kmol/kg, n counts gas only.
//! Analytic derivatives use the converged hard truncation and a fixed
//! product pool/active phase set. Finite differences use verified public solves.
//! Analytic ions, reduced element systems and condensed UV/SV states are rejected.
//! Very small finite-difference steps amplify the forward solve's constraint
//! residual, particularly when the exact constrained derivative is zero.
use crate::{
    Error, R, Result, mixture,
    solver::{Mode, Solver, State},
};
use nalgebra::DMatrix;
use std::collections::{BTreeMap, BTreeSet};
use thermo_species::Constituent;

#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct EqDerivatives {
    pub temperature: Vec<f64>,
    pub total_moles: Vec<f64>,
    pub species_amounts: BTreeMap<String, Vec<f64>>,
    pub enthalpy: Vec<f64>,
    pub energy: Vec<f64>,
    pub gibbs_energy: Vec<f64>,
    pub entropy: Vec<f64>,
    /// Maximum scaled defect in J du/dtheta + Rx (None for finite differences).
    pub closure_defect: Option<f64>,
}

fn inputs(mode: Mode) -> (f64, f64) {
    match mode {
        Mode::TP(a, b)
        | Mode::HP(a, b)
        | Mode::SP(a, b)
        | Mode::TV(a, b)
        | Mode::UV(a, b)
        | Mode::SV(a, b) => (a, b),
    }
}

fn with_inputs(mode: Mode, a: f64, b: f64) -> Mode {
    match mode {
        Mode::TP(..) => Mode::TP(a, b),
        Mode::HP(..) => Mode::HP(a, b),
        Mode::SP(..) => Mode::SP(a, b),
        Mode::TV(..) => Mode::TV(a, b),
        Mode::UV(..) => Mode::UV(a, b),
        Mode::SV(..) => Mode::SV(a, b),
    }
}

impl Solver {
    pub fn derivatives(&self, mode: Mode, weights: &[f64]) -> Result<EqDerivatives> {
        EqDerivatives::compute_derivatives(self, mode, weights)
    }

    pub fn finite_difference_derivatives(
        &self,
        mode: Mode,
        weights: &[f64],
        relative_step: f64,
    ) -> Result<EqDerivatives> {
        EqDerivatives::compute_fd(self, mode, weights, relative_step)
    }
}

impl EqDerivatives {
    /// Solve afresh so derivatives cannot refer to a stale or unrelated state.
    pub fn compute_derivatives(solver: &Solver, mode: Mode, weights: &[f64]) -> Result<Self> {
        let state = solver.solve(mode, weights)?;
        if solver.ions {
            return Err(Error::Unsupported("Analytic ion derivatives require the core's active-ion and excluded-element masks; use finite differences".into()));
        }
        let system = System::new(&state, mode, true)?;
        let sum: f64 = weights.iter().sum();
        let mut db = DMatrix::zeros(system.elements.len(), weights.len());
        for (i, e) in system.elements.iter().enumerate() {
            let b: f64 = solver
                .reactants
                .iter()
                .zip(weights)
                .map(|(r, w)| w / sum * count(&r.formula, *e) / r.molecular_weight)
                .sum();
            for (j, r) in solver.reactants.iter().enumerate() {
                db[(i, j)] = (count(&r.formula, *e) / r.molecular_weight - b) / sum;
            }
        }
        for r in &solver.reactants {
            if r.formula.iter().any(|(_, e)| !system.elements.contains(e)) {
                return Err(Error::Unsupported(
                    "Analytic weight derivatives cannot introduce an absent element".into(),
                ));
            }
        }
        system.compute(&state, mode, &db)
    }

    /// Central differences; zero weights use a second-order forward stencil.
    /// Reject product-pool/phase changes and propagate every failed solve.
    pub fn compute_fd(
        solver: &Solver,
        mode: Mode,
        weights: &[f64],
        relative_step: f64,
    ) -> Result<Self> {
        if !relative_step.is_finite() || relative_step <= 0.0 || relative_step >= 0.1 {
            return Err(Error::Input(
                "Finite-difference relative step must be between 0 and 0.1".into(),
            ));
        }
        let base = solver.solve(mode, weights)?;
        let (a, b) = inputs(mode);
        let mut out = empty(&base, weights.len() + 2);
        let sum: f64 = weights.iter().sum();
        for col in 0..weights.len() + 2 {
            let value = if col == 0 {
                a
            } else if col == 1 {
                b
            } else {
                weights[col - 2]
            };
            let scale = if col == 0 {
                match mode {
                    Mode::TP(..) | Mode::TV(..) => a.abs(),
                    Mode::HP(..) | Mode::UV(..) => a
                        .abs()
                        .max(R * base.temperature_k / base.molecular_weight_kg_kmol),
                    _ => a.abs().max(R / base.molecular_weight_kg_kmol),
                }
            } else if col == 1 {
                b
            } else if value > 0.0 {
                value
            } else {
                sum
            };
            let mut h = relative_step * scale;
            if col >= 2 && value > 0.0 {
                h = h.min(value * 0.25);
            }
            if h <= 0.0 || !h.is_finite() || value + h == value {
                return Err(Error::Input(
                    "Finite-difference step is not representable".into(),
                ));
            }
            let solve = |delta: f64| -> Result<State> {
                let mut w = weights.to_vec();
                let m = if col == 0 {
                    with_inputs(mode, a + delta, b)
                } else if col == 1 {
                    with_inputs(mode, a, b + delta)
                } else {
                    w[col - 2] += delta;
                    mode
                };
                let state = solver.solve(m, &w)?;
                if state.amounts.iter().map(|(_, s)| *s).collect::<Vec<_>>()
                    != base.amounts.iter().map(|(_, s)| *s).collect::<Vec<_>>()
                {
                    return Err(Error::Domain(
                        "Finite-difference perturbation changes the product pool or active phases"
                            .into(),
                    ));
                }
                Ok(state)
            };
            let plus = solve(h)?;
            if col >= 2 && value == 0.0 {
                let plus2 = solve(2.0 * h)?;
                difference(
                    &mut out,
                    col,
                    &[(&base, -3.0), (&plus, 4.0), (&plus2, -1.0)],
                    2.0 * h,
                );
            } else {
                let minus = solve(-h)?;
                difference(&mut out, col, &[(&plus, 1.0), (&minus, -1.0)], 2.0 * h);
            }
        }
        validate(out)
    }
}

fn count(formula: &[(f64, Constituent)], e: Constituent) -> f64 {
    formula
        .iter()
        .filter(|(_, c)| *c == e)
        .map(|(n, _)| n)
        .sum()
}

struct System {
    elements: Vec<Constituent>,
    a: DMatrix<f64>,
    j: DMatrix<f64>,
    effective: Vec<f64>,
    h: Vec<f64>,
    s: Vec<f64>,
    gas: Vec<usize>,
    condensed: Vec<usize>,
    n_col: usize,
    t_col: usize,
}

impl System {
    fn new(state: &State, mode: Mode, truncate: bool) -> Result<Self> {
        if !state.equilibrium {
            return Err(Error::Unsupported(
                "Equilibrium derivatives require an equilibrium state".into(),
            ));
        }
        let seed = state.equilibrium_seed.as_ref().ok_or_else(|| {
            Error::Unsupported("Equilibrium derivatives require the converged core state".into())
        })?;
        if seed.elements_reduced {
            return Err(Error::Unsupported(
                "Analytic derivatives after element reduction are unsupported; use finite differences".into(),
            ));
        }
        let elements: Vec<_> = state
            .amounts
            .iter()
            .filter(|(n, _)| *n > 0.0)
            .flat_map(|(_, s)| s.data().constituents().iter().map(|(_, e)| *e))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let gas: Vec<_> = state
            .amounts
            .iter()
            .enumerate()
            .filter(|(_, (_, s))| s.data().phase() == 0)
            .map(|(i, _)| i)
            .collect();
        let condensed: Vec<_> = state
            .amounts
            .iter()
            .enumerate()
            .filter(|(_, (_, s))| s.data().phase() != 0)
            .map(|(i, _)| i)
            .collect();
        if !condensed.is_empty() && matches!(mode, Mode::UV(..) | Mode::SV(..)) {
            return Err(Error::Unsupported(
                "Analytic UV/SV derivatives with condensed phases do not match public-solver finite differences; use finite differences".into(),
            ));
        }
        let ne = elements.len();
        let n_col = ne + condensed.len();
        let const_p = matches!(mode, Mode::TP(..) | Mode::HP(..) | Mode::SP(..));
        let const_t = matches!(mode, Mode::TP(..) | Mode::TV(..));
        let const_s = matches!(mode, Mode::SP(..) | Mode::SV(..));
        let t_col = n_col + usize::from(const_p);
        let size = t_col + usize::from(!const_t);
        let n = 1.0 / state.molecular_weight_kg_kmol;
        let t = state.temperature_k;
        let p = state.pressure_pa;
        let a = DMatrix::from_fn(state.amounts.len(), ne, |i, k| {
            count(state.amounts[i].1.data().constituents(), elements[k])
        });
        let effective: Vec<_> = state
            .amounts
            .iter()
            .enumerate()
            .map(|(i, (nj, s))| {
                if truncate
                    && s.data().phase() == 0
                    && seed.ln_nj[i] <= n.ln() - seed.converged_tsize
                {
                    0.0
                } else {
                    *nj
                }
            })
            .collect();
        let h: Vec<_> = state
            .amounts
            .iter()
            .map(|(_, s)| s.data().h_over_r(t) / t)
            .collect();
        let s: Vec<_> = state
            .amounts
            .iter()
            .map(|(nj, s)| {
                s.data().s_over_r(t)
                    - if s.data().phase() == 0 && *nj > 0.0 {
                        (nj * p / 1e5 / n).ln()
                    } else {
                        0.0
                    }
            })
            .collect();
        let cp: Vec<_> = state
            .amounts
            .iter()
            .map(|(_, s)| s.data().cp_over_r(t))
            .collect();
        let energy: Vec<_> = h
            .iter()
            .map(|h| h - if const_p { 0.0 } else { 1.0 })
            .collect();
        let mut j = DMatrix::zeros(size, size);
        for k in 0..ne {
            for l in 0..ne {
                j[(k, l)] = gas
                    .iter()
                    .map(|&i| effective[i] * a[(i, k)] * a[(i, l)])
                    .sum();
            }
            for (c, &i) in condensed.iter().enumerate() {
                j[(k, ne + c)] = a[(i, k)];
                j[(ne + c, k)] = a[(i, k)];
            }
            if const_p {
                j[(k, n_col)] = gas.iter().map(|&i| effective[i] * a[(i, k)]).sum();
                j[(n_col, k)] = j[(k, n_col)];
            }
            if !const_t {
                j[(k, t_col)] = gas
                    .iter()
                    .map(|&i| effective[i] * a[(i, k)] * energy[i])
                    .sum();
            }
        }
        if !const_t {
            for (c, &i) in condensed.iter().enumerate() {
                j[(ne + c, t_col)] = energy[i];
            }
        }
        if const_p {
            j[(n_col, n_col)] = gas.iter().map(|&i| effective[i]).sum::<f64>() - n;
            if !const_t {
                j[(n_col, t_col)] = gas.iter().map(|&i| effective[i] * h[i]).sum();
            }
        }
        if !const_t {
            let q: Vec<_> = if const_s { s.clone() } else { energy.clone() };
            for k in 0..ne {
                j[(t_col, k)] = gas
                    .iter()
                    .map(|&i| {
                        effective[i]
                            * a[(i, k)]
                            * (q[i] - if !const_p && const_s { 1.0 } else { 0.0 })
                    })
                    .sum();
            }
            for (c, &i) in condensed.iter().enumerate() {
                j[(t_col, ne + c)] = q[i];
            }
            if const_p {
                j[(t_col, n_col)] = gas.iter().map(|&i| effective[i] * q[i]).sum();
            }
            j[(t_col, t_col)] = effective
                .iter()
                .zip(&cp)
                .map(|(nj, cp)| nj * (cp - if const_p { 0.0 } else { 1.0 }))
                .sum::<f64>()
                + gas
                    .iter()
                    .map(|&i| {
                        effective[i]
                            * (q[i] - if !const_p && const_s { 1.0 } else { 0.0 })
                            * energy[i]
                    })
                    .sum::<f64>();
            if !const_s {
                j[(t_col, t_col)] += inputs(mode).0 / R / t
                    - effective
                        .iter()
                        .zip(&energy)
                        .map(|(nj, e)| nj * e)
                        .sum::<f64>();
            }
        }
        Ok(Self {
            elements,
            a,
            j,
            effective,
            h,
            s,
            gas,
            condensed,
            n_col,
            t_col,
        })
    }

    fn compute(&self, state: &State, mode: Mode, db: &DMatrix<f64>) -> Result<EqDerivatives> {
        let const_p = matches!(mode, Mode::TP(..) | Mode::HP(..) | Mode::SP(..));
        let const_t = matches!(mode, Mode::TP(..) | Mode::TV(..));
        let const_s = matches!(mode, Mode::SP(..) | Mode::SV(..));
        let ne = self.elements.len();
        let cols = db.ncols() + 2;
        let t = state.temperature_k;
        let n = 1.0 / state.molecular_weight_kg_kmol;
        let p = state.pressure_pa;
        let v = 1.0 / state.density_kg_m3;
        let p_fac = if const_p { 1.0 / p } else { -1.0 / v };
        let mut rhs = DMatrix::zeros(self.j.nrows(), cols);
        for k in 0..ne {
            rhs[(k, 1)] = self
                .gas
                .iter()
                .map(|&i| self.effective[i] * self.a[(i, k)])
                .sum::<f64>()
                * p_fac;
            if const_t {
                rhs[(k, 0)] = -self
                    .gas
                    .iter()
                    .map(|&i| {
                        self.effective[i]
                            * self.a[(i, k)]
                            * (self.h[i] - if const_p { 0.0 } else { 1.0 })
                    })
                    .sum::<f64>()
                    / t;
            }
            for c in 0..db.ncols() {
                rhs[(k, c + 2)] = db[(k, c)];
            }
        }
        if const_t {
            for (c, &i) in self.condensed.iter().enumerate() {
                rhs[(ne + c, 0)] = -self.h[i] / t;
            }
        }
        if const_p {
            rhs[(self.n_col, 1)] = self.gas.iter().map(|&i| self.effective[i]).sum::<f64>() / p;
            if const_t {
                rhs[(self.n_col, 0)] = -self
                    .gas
                    .iter()
                    .map(|&i| self.effective[i] * self.h[i])
                    .sum::<f64>()
                    / t;
            }
        }
        if !const_t {
            rhs[(self.t_col, 0)] = if const_s { 1.0 / R } else { 1.0 / R / t };
            rhs[(self.t_col, 1)] = self
                .gas
                .iter()
                .map(|&i| {
                    self.effective[i]
                        * (if const_s {
                            self.s[i] + if const_p { 1.0 } else { 0.0 }
                        } else {
                            self.h[i] - if const_p { 0.0 } else { 1.0 }
                        })
                })
                .sum::<f64>()
                * p_fac;
        }
        let du =
            self.j.clone().lu().solve(&rhs).ok_or_else(|| {
                Error::Convergence("Singular equilibrium derivative Jacobian".into())
            })?;
        let residual = &self.j * &du - &rhs;
        let mut defect = 0.0_f64;
        for i in 0..self.j.nrows() {
            for c in 0..cols {
                let scale = rhs[(i, c)].abs()
                    + (0..self.j.ncols())
                        .map(|k| (self.j[(i, k)] * du[(k, c)]).abs())
                        .sum::<f64>();
                defect = defect.max(residual[(i, c)].abs() / scale.max(f64::MIN_POSITIVE));
            }
        }
        if !defect.is_finite() || defect > 1e-10 {
            return Err(Error::Convergence(format!(
                "Equilibrium derivative closure defect {defect}"
            )));
        }
        let mut out = empty(state, cols);
        out.closure_defect = Some(defect);
        for c in 0..cols {
            let dt = if const_t {
                f64::from(c == 0)
            } else {
                t * du[(self.t_col, c)]
            };
            let mut dn = if const_p {
                n * du[(self.n_col, c)]
            } else {
                0.0
            };
            let dlog = if const_p {
                -dn / n + if c == 1 { 1.0 / p } else { 0.0 }
            } else {
                dt / t - if c == 1 { 1.0 / v } else { 0.0 }
            };
            let mut dh = state.cp_frozen_j_kg_k * dt;
            let mut ds = state.cp_frozen_j_kg_k / t * dt
                - R * self.gas.iter().map(|&i| state.amounts[i].0).sum::<f64>() * dlog;
            for (i, &(nj, sp)) in state.amounts.iter().enumerate() {
                let dnj = if sp.data().phase() == 0 {
                    nj * ((0..ne).map(|k| self.a[(i, k)] * du[(k, c)]).sum::<f64>()
                        + self.h[i] / t * dt
                        - dlog)
                } else {
                    let idx = self
                        .condensed
                        .iter()
                        .position(|&k| k == i)
                        .ok_or_else(|| Error::Data("Missing active condensed phase".into()))?;
                    du[(ne + idx, c)]
                };
                out.species_amounts
                    .get_mut(sp.data().symbol())
                    .ok_or_else(|| Error::Data("Missing derivative species".into()))?[c] = dnj;
                if !const_p && sp.data().phase() == 0 {
                    dn += self.effective[i] * if nj > 0.0 { dnj / nj } else { 0.0 };
                }
                dh += R * t * self.h[i] * dnj;
                ds += R * (self.s[i] - if sp.data().phase() == 0 { 1.0 } else { 0.0 }) * dnj;
            }
            if matches!(mode, Mode::HP(..)) {
                dh = f64::from(c == 0);
            }
            if const_s {
                ds = f64::from(c == 0);
            }
            let energy = if matches!(mode, Mode::UV(..)) {
                f64::from(c == 0)
            } else {
                dh - R * (n * dt + t * dn)
            };
            out.temperature[c] = dt;
            out.total_moles[c] = dn;
            out.enthalpy[c] = dh;
            out.energy[c] = energy;
            out.entropy[c] = ds;
            out.gibbs_energy[c] = dh - state.entropy_j_kg_k * dt - t * ds;
        }
        validate(out)
    }
}

fn empty(state: &State, cols: usize) -> EqDerivatives {
    EqDerivatives {
        temperature: vec![0.0; cols],
        total_moles: vec![0.0; cols],
        species_amounts: state
            .amounts
            .iter()
            .map(|(_, s)| (s.data().symbol().into(), vec![0.0; cols]))
            .collect(),
        enthalpy: vec![0.0; cols],
        energy: vec![0.0; cols],
        gibbs_energy: vec![0.0; cols],
        entropy: vec![0.0; cols],
        closure_defect: None,
    }
}

fn difference(out: &mut EqDerivatives, col: usize, states: &[(&State, f64)], denominator: f64) {
    let differentiate = |value: &dyn Fn(&State) -> f64| -> f64 {
        let base = value(states[0].0);
        states
            .iter()
            .map(|(s, c)| c * (value(s) - base))
            .sum::<f64>()
            / denominator
    };
    out.temperature[col] = differentiate(&|s| s.temperature_k);
    out.total_moles[col] = differentiate(&|s| 1.0 / s.molecular_weight_kg_kmol);
    out.enthalpy[col] = differentiate(&|s| s.enthalpy_j_kg);
    out.energy[col] = differentiate(&|s| s.energy_j_kg);
    out.gibbs_energy[col] = differentiate(&mixture::gibbs_energy);
    out.entropy[col] = differentiate(&|s| s.entropy_j_kg_k);
    for (name, row) in &mut out.species_amounts {
        row[col] = differentiate(&|s| {
            s.amounts
                .iter()
                .find(|(_, sp)| sp.data().symbol() == name)
                .map_or(0.0, |(nj, _)| *nj)
        });
    }
}

fn validate(out: EqDerivatives) -> Result<EqDerivatives> {
    if [
        &out.temperature,
        &out.total_moles,
        &out.enthalpy,
        &out.energy,
        &out.gibbs_energy,
        &out.entropy,
    ]
    .into_iter()
    .chain(out.species_amounts.values())
    .flatten()
    .any(|v| !v.is_finite())
    {
        return Err(Error::Convergence(
            "Nonfinite equilibrium derivatives".into(),
        ));
    }
    Ok(out)
}

pub(crate) fn cv_equilibrium(state: &State) -> Result<f64> {
    let mode = Mode::TP(state.temperature_k, state.pressure_pa);
    let system = System::new(state, mode, false)?;
    let d = system.compute(state, mode, &DMatrix::zeros(system.elements.len(), 0))?;
    let n = 1.0 / state.molecular_weight_kg_kmol;
    let dlnv_dlnp = state.pressure_pa / n * d.total_moles[1] - 1.0;
    let cv = -state.cp_equilibrium_j_kg_k / (state.gamma_s * dlnv_dlnp);
    if !cv.is_finite() || cv <= 0.0 {
        return Err(Error::Convergence("Nonphysical equilibrium Cv".into()));
    }
    Ok(cv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::Reactant;

    #[test]
    fn recovered_cutoff_derivatives_match_central_differences() {
        // The complete AP/Al pool at 300 K, 1 bar takes singular recovery,
        // then converges with liquid water and Al(OH)3(a), without reducing elements.
        let solver = Solver::new(
            ["NH4CLO4(I)", "AL(cr)"]
                .iter()
                .map(|name| Reactant::named(name).unwrap())
                .collect(),
            true,
            false,
            None,
            vec![],
        )
        .unwrap();
        let mode = Mode::TP(300.0, 1e5);
        let weights = [4.0, 1.0];
        let state = solver.solve(mode, &weights).unwrap();
        let seed = state.equilibrium_seed.as_ref().unwrap();
        assert_eq!(seed.converged_tsize, 80.0);
        assert!(!seed.elements_reduced);
        let analytic = solver.derivatives(mode, &weights).unwrap();
        let fd = solver
            .finite_difference_derivatives(mode, &weights, 1e-4)
            .unwrap();
        assert!(analytic.closure_defect.unwrap() < 1e-10);
        for (actual, expected) in [
            (&analytic.temperature, &fd.temperature),
            (&analytic.total_moles, &fd.total_moles),
            (&analytic.enthalpy, &fd.enthalpy),
            (&analytic.energy, &fd.energy),
            (&analytic.gibbs_energy, &fd.gibbs_energy),
            (&analytic.entropy, &fd.entropy),
        ] {
            for (col, (&a, &f)) in actual.iter().zip(expected).enumerate() {
                assert!(
                    (a - f).abs() <= 1e-5 * a.abs() + 1e-18,
                    "column {col}: analytic={a:e}, central difference={f:e}"
                );
            }
        }
        for (name, actual) in &analytic.species_amounts {
            for (col, (&a, &f)) in actual.iter().zip(&fd.species_amounts[name]).enumerate() {
                // Independent solve noise dominates derivatives of vanishing traces.
                assert!(
                    (a - f).abs() <= 1e-4 * a.abs() + 1e-12,
                    "{name} column {col}: analytic={a:e}, central difference={f:e}"
                );
            }
        }
    }
}
