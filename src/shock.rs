//! Incident and wall-reflected shocks, derived from NASA CEA `source/shock.f90`
//! at 2f79a647e85737742fb08a127da00d1ea1d069da (Apache-2.0).
//! Equilibrium states use the public solver; transport is excluded. Velocities
//! at stations 1/2 use the incident frame, and station 5 uses the reflected frame.
use crate::{
    Error, R, Result,
    solver::{Mode, Solver, State, frozen_state},
};
use std::str::FromStr;
use thermo_species::Species;

/// Specify exactly one incident speed or Mach number.
#[derive(Debug, Clone, Copy)]
pub enum ShockSpeed {
    Velocity(f64),
    Mach(f64),
}

#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object)
)]
#[derive(Debug, Clone)]
pub struct ShockStation {
    pub state: State,
    pub velocity_m_s: f64,
    pub mach: f64,
    /// NASA's reflected-frozen convention uses station 2's frozen Cp basis.
    /// The same legacy basis is used by `state.gamma_s` at that station.
    pub sonic_velocity_m_s: f64,
}

#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object)
)]
#[derive(Debug, Clone)]
pub struct ShockResult {
    pub converged: bool,
    pub last_error: i32,
    pub station1: ShockStation,
    pub station2: ShockStation,
    pub station5: Option<ShockStation>,
    pub p2_p1: f64,
    pub t2_t1: f64,
    /// Molecular-weight ratio, not a Mach ratio.
    pub m2_m1: f64,
    pub rho2_rho1: f64,
    /// Station-2 gas speed in the wall/unshocked-gas frame.
    pub v2_m_s: f64,
    pub p5_p2: Option<f64>,
    pub t5_t2: Option<f64>,
    pub m5_m2: Option<f64>,
    pub rho5_rho2: Option<f64>,
    /// Upstream gas speed in the reflected-shock frame.
    pub u5_plus_v2_m_s: Option<f64>,
}

impl ShockResult {
    fn failed(initial: State, failed: State, speed: f64, reflected: bool) -> Self {
        let mut unused = fail_state(initial.clone());
        unused.amounts.clear();
        Self {
            converged: false,
            last_error: 8,
            station1: station(initial, speed),
            station2: station(fail_state(failed), 0.0),
            station5: reflected.then(|| station(unused, 0.0)),
            p2_p1: 0.0,
            t2_t1: 0.0,
            m2_m1: 0.0,
            rho2_rho1: 0.0,
            v2_m_s: 0.0,
            p5_p2: reflected.then_some(0.0),
            t5_t2: reflected.then_some(0.0),
            m5_m2: reflected.then_some(0.0),
            rho5_rho2: reflected.then_some(0.0),
            u5_plus_v2_m_s: reflected.then_some(0.0),
        }
    }
}

fn station(state: State, velocity_m_s: f64) -> ShockStation {
    let sonic_velocity_m_s = if state.temperature_k == 0.0 {
        0.0
    } else {
        state.sound_speed()
    };
    ShockStation {
        state,
        velocity_m_s,
        mach: if sonic_velocity_m_s == 0.0 {
            0.0
        } else {
            velocity_m_s / sonic_velocity_m_s
        },
        sonic_velocity_m_s,
    }
}

// NASA clears the failed station's report but retains composition and stored M/MW/Cv_fr.
fn fail_state(mut state: State) -> State {
    if !state.equilibrium {
        // Frozen shock reports are first filled after convergence.
        state.molecular_weight_kg_kmol = 0.0;
        state.mean_molecular_weight = 0.0;
        state.cv_frozen_j_kg_k = 0.0;
        state.mole_fractions.values_mut().for_each(|v| *v = 0.0);
        state.mass_fractions.values_mut().for_each(|v| *v = 0.0);
    }
    state.converged = false;
    state.last_error = 8;
    state.post_processed = false;
    state.temperature_k = 0.0;
    state.gas_moles = 0.0;
    state.pressure_pa = 0.0;
    state.density_kg_m3 = 0.0;
    state.enthalpy_j_kg = 0.0;
    state.energy_j_kg = 0.0;
    state.gibbs_energy_j_kg = 0.0;
    state.entropy_j_kg_k = 0.0;
    state.cp_frozen_j_kg_k = 0.0;
    state.cp_equilibrium_j_kg_k = 0.0;
    state.cv_equilibrium_j_kg_k = 0.0;
    state.gamma_s = 0.0;
    state.dlnv_dlnp = 0.0;
    state.dlnv_dlnt = 0.0;
    state
}

/// NASA's analytic log-pressure/log-temperature shock equations.
#[allow(clippy::too_many_arguments)]
fn solve_jump(
    up: &State,
    speed: f64,
    reflected: bool,
    frozen: bool,
    mut p_ratio: f64,
    mut t_ratio: f64,
    mut state_at: impl FnMut(f64, f64) -> Result<State>,
) -> Result<(State, i32, f64, f64, f64)> {
    let wm = up.molecular_weight_kg_kmol;
    let t0 = up.temperature_k;
    let mu = wm * speed * speed / (R * t0);
    let h_init = up.enthalpy_j_kg / R;
    let mut last = None;
    let mut last_valid = None;
    for iter in 1..=60 {
        let t = t_ratio * t0;
        let st = state_at(t, up.pressure_pa * p_ratio)?;
        if !frozen
            && (!(st.temperature_k > 0.0 && st.gas_moles > 0.0) || (reflected && !st.converged))
        {
            return Ok(if !reflected && let Some(last_valid) = last_valid {
                last_valid
            } else {
                (fail_state(st), 8, 0.0, 0.0, 0.0)
            });
        }
        let wmk = 1.0 / st.gas_moles;
        let rho12 = wm * t_ratio / (p_ratio * wmk);
        last_valid = Some((st.clone(), 9, p_ratio, t_ratio, rho12));
        let cp = st.cp_equilibrium_j_kg_k / R;
        let (vp, vt) = (st.dlnv_dlnp, st.dlnv_dlnt);
        let rho52 = 1.0 / rho12;
        let mut tmp = if reflected {
            -mu * rho52 / (rho52 - 1.0).powi(2)
        } else {
            rho12 * mu
        };
        let a = (-tmp).mul_add(vp, -p_ratio);
        let b = -tmp * vt;
        let f = if reflected {
            p_ratio - 1.0 + tmp * (rho52 - 1.0)
        } else {
            p_ratio - 1.0 + tmp - mu
        };
        tmp *= t0 / wm;
        if !reflected {
            tmp *= rho12;
        }
        let c = (-tmp).mul_add(vp, t * (vt - 1.0) / wmk);
        let d = (-tmp).mul_add(vt, -t * cp);
        let g = st.iteration_enthalpy() / R
            - h_init
            - speed
                * speed
                * if reflected {
                    (rho52 + 1.0) / (rho52 - 1.0)
                } else {
                    (-rho12).mul_add(rho12, 1.0)
                }
                / (2.0 * R);
        // NASA's fused products retain the determinant on singular ion-shock paths.
        let det = a.mul_add(d, -b * c);
        let mut dp = f.mul_add(d, -g * b) / det;
        let mut dt = a.mul_add(g, -c * f) / det;
        let correction = dp.abs().max(dt.abs());
        if correction < 5e-5 {
            return Ok((st, 0, p_ratio, t_ratio, rho12));
        }
        let cap = if iter > 20 {
            0.04879016
        } else if iter > 12 {
            0.09531018
        } else if iter > 4 {
            0.22314355
        } else {
            0.40546511
        };
        if correction / cap > 1.0 {
            dp /= correction / cap;
            dt /= correction / cap;
        }
        p_ratio = (p_ratio.ln() + dp).exp();
        t_ratio = (t_ratio.ln() + dt).exp();
        if iter == 1 && t_ratio >= 21000.0 / t0 {
            return Ok((fail_state(st), 8, 0.0, 0.0, 0.0));
        }
        last = Some((st, 9, p_ratio, t_ratio, rho12));
    }
    let last = last.unwrap();
    Ok(if reflected || frozen {
        (fail_state(last.0), 8, 0.0, 0.0, 0.0)
    } else {
        last
    })
}

impl Solver {
    /// Solve a gaseous incident shock and optionally its wall reflection.
    /// Weights are masses; temperature K, pressure Pa, velocity m/s.
    /// All four incident/reflected equilibrium/frozen combinations are allowed.
    #[allow(clippy::too_many_arguments)]
    pub fn shock(
        &self,
        weights: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        speed: ShockSpeed,
        reflected: bool,
        incident_frozen: bool,
        reflected_frozen: bool,
    ) -> Result<ShockResult> {
        self.feeds(weights)?;
        let total: f64 = weights.iter().sum();
        let mut amounts = vec![];
        for (r, w) in self
            .reactants
            .iter()
            .zip(weights)
            .filter(|(_, w)| **w > 0.0)
        {
            let sp = Species::from_str(&r.name)
                .map_err(|_| Error::Input("Shock reactants must be NASA gas species".into()))?;
            if sp.data().phase() != 0
                || r.properties(temperature_k)?.0.is_none()
                || r.molecular_weight != sp.data().mw()
                || r.formula != sp.data().constituents()
            {
                return Err(Error::Input(
                    "Shock reactants must be unmodified NASA gas species with heat capacity".into(),
                ));
            }
            if let Some((n, _)) = amounts.iter_mut().find(|(_, s)| *s == sp) {
                *n += w / total / r.molecular_weight;
            } else {
                amounts.push((w / total / r.molecular_weight, sp));
            }
        }
        if !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !pressure_pa.is_finite()
            || pressure_pa <= 0.0
        {
            return Err(Error::Domain(
                "Shock temperature/pressure outside the thermodynamic domain".into(),
            ));
        }
        let initial = frozen_state(&amounts, temperature_k, pressure_pa, None)?;
        let u1 = match speed {
            ShockSpeed::Velocity(u) => u,
            ShockSpeed::Mach(m) => m * initial.sound_speed(),
        };
        let mach1 = u1 / initial.sound_speed();
        if !mach1.is_finite() || mach1 <= 1.0 {
            return Err(Error::Input(
                "Incident shock must have finite Mach number > 1".into(),
            ));
        }
        let gamma = initial.gamma_s;
        let p21 = (2.0 * gamma * mach1 * mach1 - gamma + 1.0) / (gamma + 1.0);
        let mut seed = self.initial_seed(weights, temperature_k, Some(&initial.amounts))?;
        let t21 = if incident_frozen {
            p21 * (2.0 / (mach1 * mach1) + gamma - 1.0) / (gamma + 1.0)
        } else {
            let hp = self.solve_with_start(
                Mode::HP(initial.enthalpy_j_kg + u1 * u1 / 2.0, pressure_pa * p21),
                weights,
                Some(&seed),
            )?;
            if hp.converged {
                seed = hp.equilibrium_seed.clone().unwrap();
            }
            if !(hp.temperature_k > 0.0 && hp.gas_moles > 0.0) {
                return Ok(ShockResult::failed(initial, hp, u1, reflected));
            }
            hp.temperature_k / temperature_k
        }
        .min(21000.0 / temperature_k);
        let (mut incident, status, p21, t21, rho12) =
            solve_jump(&initial, u1, false, incident_frozen, p21, t21, |t, p| {
                if incident_frozen {
                    initial.frozen(t, p)
                } else {
                    let st = self.solve_state(Mode::TP(t, p), weights, Some(&seed), true)?;
                    if st.converged {
                        seed = st.equilibrium_seed.clone().unwrap();
                    }
                    Ok(st)
                }
            })?;
        if status == 8 {
            return Ok(ShockResult::failed(initial, incident, u1, reflected));
        }
        incident.pressure_pa = p21 * pressure_pa;
        incident.transport_atomic =
            !incident.converged || (!incident_frozen && reflected && reflected_frozen);
        let u2 = u1 * rho12;
        let v2 = u1 - u2;
        let mut result = ShockResult {
            converged: status == 0,
            last_error: status,
            p2_p1: p21,
            t2_t1: t21,
            m2_m1: incident.molecular_weight_kg_kmol / initial.molecular_weight_kg_kmol,
            rho2_rho1: 1.0 / rho12,
            v2_m_s: v2,
            station1: station(initial, u1),
            station2: station(incident, u2),
            station5: None,
            p5_p2: None,
            t5_t2: None,
            m5_m2: None,
            rho5_rho2: None,
            u5_plus_v2_m_s: None,
        };
        if reflected && result.station2.state.temperature_k > 0.0 {
            let up = &result.station2.state;
            let mu = up.density_kg_m3 * v2 * v2 / up.pressure_pa;
            let b = (1.0 + mu + 2.0) / 2.0;
            let p52 = b + (b * b - 2.0).sqrt();
            let mut seed = self.initial_seed(weights, up.temperature_k, Some(&up.amounts))?;
            let (mut down, status, p52, t52, rho25) = solve_jump(
                up,
                v2,
                true,
                reflected_frozen,
                p52,
                2.0_f64.min(21000.0 / up.temperature_k),
                |t, p| {
                    if reflected_frozen {
                        up.frozen(t, p)
                    } else {
                        let st = self.solve_state(Mode::TP(t, p), weights, Some(&seed), false)?;
                        if st.converged {
                            seed = st.equilibrium_seed.clone().unwrap();
                        }
                        Ok(st)
                    }
                },
            )?;
            result.converged = status == 0;
            result.last_error = status;
            if status == 8 {
                result.p5_p2 = Some(0.0);
                result.t5_t2 = Some(0.0);
                result.m5_m2 = Some(0.0);
                result.rho5_rho2 = Some(0.0);
                result.u5_plus_v2_m_s = Some(0.0);
                result.station5 = Some(station(down, 0.0));
                return Ok(result);
            }
            let ratio = 1.0 / rho25;
            let u5 = v2 / (ratio - 1.0);
            if reflected_frozen {
                // NASA's reflected-frozen header retains incident Cp, even at T5.
                let cp = up.cp_frozen_j_kg_k;
                down.gamma_s = cp / (cp - R / up.molecular_weight_kg_kmol);
                down.cp_equilibrium_j_kg_k = cp;
                down.cv_equilibrium_j_kg_k = cp - down.gas_moles * R;
                down.cp_frozen_j_kg_k = 0.0;
                down.cv_frozen_j_kg_k = 0.0;
            }
            result.p5_p2 = Some(p52);
            result.t5_t2 = Some(t52);
            result.m5_m2 = Some(down.molecular_weight_kg_kmol / up.molecular_weight_kg_kmol);
            result.rho5_rho2 = Some(ratio);
            result.u5_plus_v2_m_s = Some(u5 + v2);
            result.station5 = Some(station(down, u5));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::Reactant;

    #[test]
    fn failed_ion_shock_partials_match_native() {
        // NASA's first TP iterate after the failed N2/CH4 shock HP initialization.
        let mut solver = Solver::new(
            ["N2", "CH4"]
                .map(|name| Reactant::named(name).unwrap())
                .to_vec(),
            true,
            true,
            None,
            vec![],
        )
        .unwrap()
        .with_trace(1e-10)
        .unwrap();
        solver.transport = true;
        let weights = solver.moles_to_weights(&[0.985, 0.015]).unwrap();
        let total: f64 = weights.iter().sum();
        let amounts: Vec<_> = solver
            .reactants
            .iter()
            .zip(&weights)
            .map(|(r, w)| {
                (
                    w / total / r.molecular_weight,
                    Species::from_str(&r.name).unwrap(),
                )
            })
            .collect();
        let seed = solver
            .initial_seed(&weights, 291.0, Some(&amounts))
            .unwrap();
        let state = solver
            .solve_state(
                Mode::TP(298.14999999999895, 407.689_193_099_219),
                &weights,
                Some(&seed),
                true,
            )
            .unwrap();
        assert!(!state.converged);
        for (actual, expected) in [
            (state.gas_moles, 19742862968.172836),
            (state.dlnv_dlnp, -530396717550.6311),
            (state.dlnv_dlnt, 0.9953529048453027),
            (state.cp_equilibrium_j_kg_k, 410380578944840.75),
        ] {
            assert!(
                (actual - expected).abs() / expected.abs() < 1e-4,
                "{actual} != {expected}"
            );
        }

        // Reproduce the first five TP solves; the unfused path produces a nonfinite sixth input.
        let initial = frozen_state(&amounts, 291.0, 29.864211683161557, None).unwrap();
        let gamma = initial.gamma_s;
        let mach = 1200.0 / initial.sound_speed();
        let p21 = (2.0 * gamma * mach * mach - gamma + 1.0) / (gamma + 1.0);
        let mut calls = 0;
        let (_, status, _, _, _) = solve_jump(
            &initial,
            1200.0,
            false,
            false,
            p21,
            298.14999999999895 / 291.0,
            |t, p| {
                calls += 1;
                assert!(t.is_finite() && t > 0.0 && p.is_finite() && p > 0.0);
                if calls <= 5 {
                    solver.solve_state(Mode::TP(t, p), &weights, Some(&seed), true)
                } else {
                    Ok(fail_state(initial.clone()))
                }
            },
        )
        .unwrap();
        assert_eq!(calls, 6);
        assert_eq!(status, 9);
    }
}
