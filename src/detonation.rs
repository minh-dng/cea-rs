//! Chapman-Jouguet equilibrium detonation, with SI inputs and outputs.
//!
//! Derived from NASA CEA's `source/detonation.f90` at commit
//! 2f79a647e85737742fb08a127da00d1ea1d069da, under Apache-2.0.
//! The initial gas retains its reactant composition. Burned states use the
//! public equilibrium solver; transport and frozen detonation are unsupported.
use crate::{
    Error, R, Result,
    solver::{Mode, Solver, State},
};
use std::str::FromStr;
use thermo_species::Species;

#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object)
)]
#[derive(Debug, Clone)]
pub struct UnburnedState {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub density_kg_m3: f64,
    pub enthalpy_j_kg: f64,
    pub molecular_weight_kg_kmol: f64,
    pub gamma_s: f64,
    pub sound_speed_m_s: f64,
}

#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object)
)]
#[derive(Debug, Clone)]
pub struct DetonationResult {
    pub converged: bool,
    pub last_error: i32,
    pub sonic_velocity_m_s: f64,
    pub unburned: UnburnedState,
    pub burned: State,
    pub pressure_ratio: f64,
    pub temperature_ratio: f64,
    pub molecular_weight_ratio: f64,
    pub density_ratio: f64,
    /// Detonation velocity divided by the unburned sound speed.
    pub mach: f64,
    /// Wave velocity relative to the stationary unburned gas, m/s.
    pub velocity_m_s: f64,
}

impl Solver {
    /// Solve an equilibrium CJ wave for mass weights and unburned K, Pa.
    ///
    /// Reactants must be named gas species with NASA heat-capacity fits.
    /// As in NASA CEA 3.3.4, `frozen=true` is unsupported; it does not mean
    /// freezing the unburned gas, whose composition is always held fixed.
    pub fn detonation(
        &self,
        weights: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        frozen: bool,
    ) -> Result<DetonationResult> {
        if frozen {
            return Err(Error::Unsupported(
                "Frozen detonation is not implemented in NASA CEA 3.3.4".into(),
            ));
        }
        self.feeds(weights)?;
        if !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !pressure_pa.is_finite()
            || pressure_pa <= 0.0
        {
            return Err(Error::Domain(
                "Unburned temperature and pressure must be positive and finite".into(),
            ));
        }
        let total: f64 = weights.iter().sum();
        let mut moles = 0.0;
        let mut enthalpy = 0.0;
        let mut cp1 = 0.0;
        for (reactant, &weight) in self.reactants.iter().zip(weights) {
            if weight == 0.0 {
                continue;
            }
            let phase = Species::from_str(&reactant.name)
                .ok()
                .map(|s| s.data().phase())
                .or_else(|| {
                    // Reactant-only fits (including Air) are absent from Species.
                    let mut records = include_str!("../vendor/thermo-species/data/thermo.inp")
                        .lines()
                        .skip_while(|line| *line != "END PRODUCTS")
                        .skip(1);
                    records.find(|line| line.get(..18).map(str::trim) == Some(&reactant.name))?;
                    records.next()?.get(51..52)?.trim().parse().ok()
                });
            if phase != Some(0) {
                return Err(Error::Input(format!(
                    "Detonation reactant requires a named gas species fit: {}",
                    reactant.name
                )));
            }
            let (cp, h, _) = reactant.properties(temperature_k)?;
            let amount = weight / total / reactant.molecular_weight;
            moles += amount;
            enthalpy += amount * h * 1000.0;
            cp1 += amount
                * cp.ok_or_else(|| {
                    Error::Unsupported("Detonation requires reactant heat capacities".into())
                })?
                * 1000.0;
        }
        let mw1 = 1.0 / moles;
        let gamma1 = cp1 / (cp1 - R / mw1);
        let rho1 = pressure_pa * mw1 / (R * temperature_k);
        let sound1 = (gamma1 * R * temperature_k / mw1).sqrt();
        if [mw1, gamma1, rho1, sound1]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || gamma1 <= 1.0
            || !enthalpy.is_finite()
        {
            return Err(Error::Domain("Nonphysical unburned gas properties".into()));
        }
        let unburned = UnburnedState {
            temperature_k,
            pressure_pa,
            density_kg_m3: rho1,
            enthalpy_j_kg: enthalpy,
            molecular_weight_kg_kmol: mw1,
            gamma_s: gamma1,
            sound_speed_m_s: sound1,
        };

        // NASA's HP estimate and three recursion steps initialize its Newton solve.
        let mut pp1 = 15.0;
        let initial = self.solve(
            Mode::HP(
                enthalpy + R * 0.75 * temperature_k * pp1 / mw1,
                pressure_pa * pp1,
            ),
            weights,
        )?;
        let gamma = initial.gamma_s;
        let cp = initial.cp_equilibrium_j_kg_k / R;
        let mut tt1 = initial.temperature_k / temperature_k;
        let temp = tt1 - 0.75 * pp1 / (cp * mw1);
        let amm = initial.molecular_weight_kg_kmol / mw1;
        for _ in 0..3 {
            let alfa = amm / tt1;
            pp1 = (1.0 + gamma) * (1.0 + (1.0 - 4.0 * gamma * alfa / (1.0 + gamma).powi(2)).sqrt())
                / (2.0 * gamma * alfa);
            let rk = pp1 * alfa;
            tt1 = temp + 0.5 * pp1 * gamma * (rk * rk - 1.0) / (mw1 * cp * rk);
        }

        let mut seed = initial.equilibrium_seed.clone();
        let mut result = None;
        for _ in 0..8 {
            let t2 = temperature_k * tt1;
            let p2 = pressure_pa * pp1;
            let mut state = self.solve_with_start(Mode::TP(t2, p2), weights, seed.as_ref())?;
            if state.converged {
                seed.clone_from(&state.equilibrium_seed);
            }
            let gamma = state.gamma_s;
            let cp = state.cp_equilibrium_j_kg_k / R;
            let mw = 1.0 / state.gas_moles;
            let rr1 = pp1 * mw / mw1 / tt1;
            let (vt, vp) = (state.dlnv_dlnt, state.dlnv_dlnp);
            let a11 = 1.0 / pp1 + gamma * rr1 * vp;
            let a12 = gamma * rr1 * vt;
            let a21 = 0.5 * gamma * (rr1 * rr1 - 1.0 - vp * (1.0 + rr1 * rr1)) + vt - 1.0;
            let a22 = -0.5 * gamma * vt * (rr1 * rr1 + 1.0) - mw * cp;
            let b1 = 1.0 / pp1 - 1.0 + gamma * (rr1 - 1.0);
            let b2 = mw * (state.iteration_enthalpy() - enthalpy) / (R * t2)
                - 0.5 * gamma * (rr1 * rr1 - 1.0);
            let determinant = a11 * a22 - a12 * a21;
            let x1 = (a22 * b1 - a12 * b2) / determinant;
            let x2 = (a11 * b2 - a21 * b1) / determinant;
            let change = x1.abs().max(x2.abs());
            let damping = if change > 0.4054652 {
                0.4054652 / change
            } else {
                1.0
            };
            pp1 *= (x1 * damping).exp();
            tt1 *= (x2 * damping).exp();
            // NASA updates T and wave speeds while retaining the preceding thermodynamic fields.
            state.temperature_k = temperature_k * tt1;
            let sonic = (gamma * R * state.temperature_k / mw).sqrt();
            let velocity = rr1 * sonic;
            let converged = change < 0.5e-4;
            result = Some(DetonationResult {
                converged,
                last_error: if converged { 0 } else { 8 },
                sonic_velocity_m_s: sonic,
                pressure_ratio: pp1,
                temperature_ratio: tt1,
                molecular_weight_ratio: mw / mw1,
                density_ratio: rr1,
                mach: velocity / unburned.sound_speed_m_s,
                velocity_m_s: velocity,
                unburned: unburned.clone(),
                burned: state,
            });
            if converged {
                break;
            }
        }
        Ok(result.unwrap())
    }
}
