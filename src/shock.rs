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

fn station(state: State, velocity_m_s: f64) -> ShockStation {
    let sonic_velocity_m_s = state.sound_speed();
    ShockStation {
        state,
        velocity_m_s,
        mach: velocity_m_s / sonic_velocity_m_s,
        sonic_velocity_m_s,
    }
}

/// NASA's analytic log-pressure/log-temperature shock equations.
fn solve_jump(
    up: &State,
    speed: f64,
    reflected: bool,
    mut p_ratio: f64,
    mut t_ratio: f64,
    mut state_at: impl FnMut(f64, f64) -> Result<State>,
) -> Result<(State, bool, f64, f64, f64)> {
    let wm = up.molecular_weight_kg_kmol;
    let t0 = up.temperature_k;
    let mu = wm * speed * speed / (R * t0);
    let h_init = up.enthalpy_j_kg / R;
    let mut last = None;
    for iter in 1..=60 {
        let t = t_ratio * t0;
        let st = state_at(t, up.pressure_pa * p_ratio)?;
        let wmk = 1.0 / st.gas_moles;
        let rho12 = wm * t_ratio / (p_ratio * wmk);
        let cp = st.cp_equilibrium_j_kg_k / R;
        let (vp, vt) = (st.dlnv_dlnp, st.dlnv_dlnt);
        let rho52 = 1.0 / rho12;
        let mut tmp = if reflected {
            -mu * rho52 / (rho52 - 1.0).powi(2)
        } else {
            rho12 * mu
        };
        let a = -tmp * vp - p_ratio;
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
        let c = -tmp * vp + t * (vt - 1.0) / wmk;
        let d = -tmp * vt - t * cp;
        let g = st.iteration_enthalpy() / R
            - h_init
            - speed * speed / (2.0 * R)
                * if reflected {
                    (rho52 + 1.0) / (rho52 - 1.0)
                } else {
                    1.0 - rho12 * rho12
                };
        let det = a * d - b * c;
        let mut dp = (f * d - g * b) / det;
        let mut dt = (a * g - c * f) / det;
        let correction = dp.abs().max(dt.abs());
        if correction < 5e-5 {
            return Ok((st, true, p_ratio, t_ratio, rho12));
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
        last = Some((st, false, p_ratio, t_ratio, rho12));
    }
    Ok(last.unwrap())
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
            hp.temperature_k / temperature_k
        }
        .min(21000.0 / temperature_k);
        let (mut incident, converged, p21, t21, rho12) =
            solve_jump(&initial, u1, false, p21, t21, |t, p| {
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
        incident.pressure_pa = p21 * pressure_pa;
        incident.transport_atomic = !incident.converged;
        let u2 = u1 * rho12;
        let v2 = u1 - u2;
        let mut result = ShockResult {
            converged,
            last_error: if converged { 0 } else { 9 },
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
            let (mut down, converged, p52, t52, rho25) =
                solve_jump(up, v2, true, p52, 2.0, |t, p| {
                    if reflected_frozen {
                        up.frozen(t, p)
                    } else {
                        let st = self.solve_with_start(Mode::TP(t, p), weights, Some(&seed))?;
                        if st.converged {
                            seed = st.equilibrium_seed.clone().unwrap();
                        }
                        Ok(st)
                    }
                })?;
            result.converged = converged;
            result.last_error = if converged { 0 } else { 8 };
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
