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

fn verify_jump(up: &State, down: &State, u: f64, v: f64) -> Result<()> {
    let mass = up.density_kg_m3 * u;
    let momentum = up.pressure_pa + mass * u;
    let energy = up.enthalpy_j_kg + u * u / 2.0;
    let errors = [
        (mass - down.density_kg_m3 * v).abs() / mass.abs().max(1.0),
        (momentum - down.pressure_pa - down.density_kg_m3 * v * v).abs() / momentum.abs().max(1.0),
        (energy - down.enthalpy_j_kg - v * v / 2.0).abs()
            / energy
                .abs()
                .max(u * u)
                .max(R * up.temperature_k / up.molecular_weight_kg_kmol)
                .max(1.0),
    ];
    if errors.iter().any(|e| !e.is_finite() || *e > 5e-5)
        || down.pressure_pa <= up.pressure_pa
        || down.density_kg_m3 <= up.density_kg_m3
        || !u.is_finite()
        || !v.is_finite()
        || v <= 0.0
        || u <= v
    {
        return Err(Error::Convergence(format!(
            "Shock failed mass, momentum or energy closure: {errors:?}"
        )));
    }
    Ok(())
}

/// NASA's two log-pressure/log-temperature Newton equations and damping.
/// Numerical TP derivatives keep the equilibrium core behind its public API.
fn solve_jump(
    up: &State,
    speed: f64,
    reflected: bool,
    mut p_ratio: f64,
    mut t_ratio: f64,
    state_at: impl Fn(f64, f64) -> Result<State>,
) -> Result<State> {
    let residual = |st: &State| {
        let d = st.density_kg_m3 / up.density_kg_m3;
        let (dynamic, kinetic) = if reflected {
            (
                up.density_kg_m3 * speed * speed * d / (d - 1.0),
                speed * speed * (d + 1.0) / (d - 1.0) / 2.0,
            )
        } else {
            (
                up.density_kg_m3 * speed * speed * (1.0 - 1.0 / d),
                speed * speed * (1.0 - 1.0 / (d * d)) / 2.0,
            )
        };
        [
            (st.pressure_pa - up.pressure_pa - dynamic) / up.pressure_pa,
            (st.enthalpy_j_kg - up.enthalpy_j_kg - kinetic) / R,
        ]
    };
    for iter in 1..=60 {
        let p = up.pressure_pa * p_ratio;
        let t = up.temperature_k * t_ratio;
        let mut st = state_at(t, p)?;
        let f = residual(&st);
        let step: f64 = 1e-4;
        let fp = residual(&state_at(t, p * step.exp())?);
        let ft = residual(&state_at(t * step.exp(), p)?);
        let a = (fp[0] - f[0]) / step;
        let b = (ft[0] - f[0]) / step;
        let c = (fp[1] - f[1]) / step;
        let d = (ft[1] - f[1]) / step;
        let det = a * d - b * c;
        let mut dp = (b * f[1] - d * f[0]) / det;
        let mut dt = (c * f[0] - a * f[1]) / det;
        if !dp.is_finite() || !dt.is_finite() {
            return Err(Error::Convergence("Singular shock Newton system".into()));
        }
        let correction = dp.abs().max(dt.abs());
        if correction < 5e-5 {
            let ratio = st.density_kg_m3 / up.density_kg_m3;
            let (u, v) = if reflected {
                (speed * ratio / (ratio - 1.0), speed / (ratio - 1.0))
            } else {
                (speed, speed / ratio)
            };
            verify_jump(up, &st, u, v)?;
            st.iterations += iter;
            return Ok(st);
        }
        let cap = if iter > 20 {
            1.05_f64.ln()
        } else if iter > 12 {
            1.1_f64.ln()
        } else if iter > 4 {
            1.25_f64.ln()
        } else {
            1.5_f64.ln()
        };
        let damping = (cap / correction).min(1.0);
        dp *= damping;
        dt *= damping;
        p_ratio *= dp.exp();
        t_ratio *= dt.exp();
    }
    Err(Error::Convergence("Shock exceeded 60 iterations".into()))
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
        let initial = frozen_state(&amounts, temperature_k, pressure_pa)?;
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
        let t21 = if incident_frozen {
            p21 * (2.0 / (mach1 * mach1) + gamma - 1.0) / (gamma + 1.0)
        } else {
            self.solve(
                Mode::HP(initial.enthalpy_j_kg + u1 * u1 / 2.0, pressure_pa * p21),
                weights,
            )?
            .temperature_k
                / temperature_k
        };
        let incident = solve_jump(&initial, u1, false, p21, t21, |t, p| {
            if incident_frozen {
                frozen_state(&amounts, t, p)
            } else {
                self.solve(Mode::TP(t, p), weights)
            }
        })?;
        let u2 = u1 * initial.density_kg_m3 / incident.density_kg_m3;
        let v2 = u1 - u2;
        let mut result = ShockResult {
            p2_p1: incident.pressure_pa / pressure_pa,
            t2_t1: incident.temperature_k / temperature_k,
            m2_m1: incident.molecular_weight_kg_kmol / initial.molecular_weight_kg_kmol,
            rho2_rho1: incident.density_kg_m3 / initial.density_kg_m3,
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
        if reflected {
            let up = &result.station2.state;
            let mu = up.density_kg_m3 * v2 * v2 / up.pressure_pa;
            let b = (1.0 + mu + 2.0) / 2.0;
            let p52 = b + (b * b - 2.0).sqrt();
            let mut down = solve_jump(up, v2, true, p52, 2.0, |t, p| {
                if reflected_frozen {
                    let mut st = frozen_state(&up.amounts, t, p)?;
                    st.element_residual = up.element_residual;
                    Ok(st)
                } else {
                    self.solve(Mode::TP(t, p), weights)
                }
            })?;
            let ratio = down.density_kg_m3 / up.density_kg_m3;
            let u5 = v2 / (ratio - 1.0);
            if reflected_frozen {
                // NASA's reflected-frozen header retains incident Cp, even at T5.
                let cp = up.cp_frozen_j_kg_k;
                down.gamma_s = cp / (cp - R / up.molecular_weight_kg_kmol);
            }
            result.p5_p2 = Some(down.pressure_pa / up.pressure_pa);
            result.t5_t2 = Some(down.temperature_k / up.temperature_k);
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
    fn conservation_guard_rejects_corrupted_mass_momentum_and_energy() {
        let solver = Solver::new(
            vec![Reactant::named("N2").unwrap()],
            false,
            false,
            None,
            vec![],
        )
        .unwrap();
        let result = solver
            .shock(
                &[1.0],
                300.0,
                1e5,
                ShockSpeed::Mach(3.0),
                false,
                true,
                false,
            )
            .unwrap();
        let up = &result.station1;
        let down = &result.station2;
        for property in 0..3 {
            let mut bad = down.state.clone();
            match property {
                0 => bad.density_kg_m3 *= 1.1,
                1 => bad.pressure_pa *= 1.1,
                _ => bad.enthalpy_j_kg += 1e5,
            }
            assert!(matches!(
                verify_jump(&up.state, &bad, up.velocity_m_s, down.velocity_m_s),
                Err(Error::Convergence(_))
            ));
        }
    }
}
