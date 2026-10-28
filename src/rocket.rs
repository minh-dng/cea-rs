//! Infinite-area chamber and isentropic nozzle expansion in native Rust.
use crate::{
    Error, R, Result,
    solver::{Mode, Solver, State},
    thermo::check_species_temperature,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowMode {
    Equilibrium,
    FrozenAtChamber,
    FrozenAtThroat,
}
#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct Station {
    pub state: State,
    pub velocity_m_s: f64,
    pub mach: f64,
    pub area_ratio: Option<f64>,
}
#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct RocketResult {
    pub chamber: Station,
    pub throat: Station,
    pub exit: Station,
    pub cstar_m_s: f64,
    pub vacuum_velocity_m_s: f64,
}

fn bracketed_root(
    mut lo: f64,
    mut hi: f64,
    mut value: impl FnMut(f64) -> Result<f64>,
) -> Result<f64> {
    let mut flo = value(lo)?;
    let mut fhi = value(hi)?;
    if !flo.is_finite() || !fhi.is_finite() {
        return Err(Error::Convergence(
            "Nonfinite scalar bracket residual".into(),
        ));
    }
    if flo == 0.0 {
        return Ok(lo);
    }
    if fhi == 0.0 {
        return Ok(hi);
    }
    if flo.signum() == fhi.signum() {
        return Err(Error::Convergence("Scalar root is not bracketed".into()));
    }
    let mut previous_side = 0;
    for _ in 0..80 {
        let candidate = lo - flo / (fhi - flo) * (hi - lo);
        let mid = if candidate.is_finite() && candidate > lo && candidate < hi {
            candidate
        } else {
            lo + (hi - lo) / 2.0
        };
        let fm = value(mid)?;
        if !fm.is_finite() {
            return Err(Error::Convergence("Nonfinite scalar residual".into()));
        }
        if fm.abs() < 1e-9 || (hi - lo).abs() < 1e-10 {
            return Ok(mid);
        }
        if fm.signum() == flo.signum() {
            lo = mid;
            flo = fm;
            // Illinois damping prevents a stationary endpoint from stalling.
            if previous_side == 1 {
                fhi *= 0.5;
            }
            previous_side = 1;
        } else {
            hi = mid;
            fhi = fm;
            if previous_side == -1 {
                flo *= 0.5;
            }
            previous_side = -1;
        }
    }
    Err(Error::Convergence(
        "Scalar root exceeded 80 iterations".into(),
    ))
}

fn frozen_state(reference: &State, t: f64, p: f64) -> Result<State> {
    let mut out = reference.clone();
    let ng: f64 = reference
        .amounts
        .iter()
        .filter(|(_, s)| s.data().phase() == 0)
        .map(|(n, _)| n)
        .sum();
    let total: f64 = reference.amounts.iter().map(|(n, _)| n).sum();
    out.trace_extrapolated_species.clear();
    let mut h = 0.0;
    let mut s = 0.0;
    let mut cp = 0.0;
    for &(n, sp) in &reference.amounts {
        if n == 0.0 {
            continue;
        }
        if n / total > 1e-11 {
            check_species_temperature(sp, t)?;
        } else if check_species_temperature(sp, t).is_err() {
            out.trace_extrapolated_species
                .push(sp.data().symbol().into());
        }
        let d = sp.data();
        h += n * d.h_over_r(t) * R;
        cp += n * d.cp_over_r(t) * R;
        let mixing = if d.phase() == 0 {
            (n / ng * p / 1e5).ln()
        } else {
            0.0
        };
        s += n * (d.s_over_r(t) - mixing) * R;
    }
    let cv = cp - ng * R;
    if !h.is_finite() || !s.is_finite() || !cv.is_finite() || cv <= 0.0 || ng <= 0.0 {
        return Err(Error::Convergence("Nonphysical frozen state".into()));
    }
    out.temperature_k = t;
    out.pressure_pa = p;
    out.enthalpy_j_kg = h;
    out.entropy_j_kg_k = s;
    out.energy_j_kg = h - ng * R * t;
    out.cp_frozen_j_kg_k = cp;
    out.cp_equilibrium_j_kg_k = cp;
    out.gamma_s = cp / cv;
    out.density_kg_m3 = p / (ng * R * t);
    out.molecular_weight_kg_kmol = 1.0 / ng;
    out.iterations = 0;
    out.constraint_residual = 0.0;
    out.equilibrium = false;
    Ok(out)
}

fn frozen_expansion(reference: &State, p: f64) -> Result<State> {
    // No phase changes or extrapolation are allowed for appreciable frozen species.
    let total: f64 = reference.amounts.iter().map(|(n, _)| n).sum();
    let min_t = reference
        .amounts
        .iter()
        .filter(|(n, _)| n / total > 1e-11)
        .map(|(_, sp)| sp.data().valid_temperature_range().0)
        .fold(300.0_f64, f64::max);
    let target = reference.entropy_j_kg_k;
    let t = bracketed_root(min_t.ln(), reference.temperature_k.ln(), |x| {
        Ok((frozen_state(reference, x.exp(), p)?.entropy_j_kg_k - target) / target.abs().max(1.0))
    })?
    .exp();
    let mut state = frozen_state(reference, t, p)?;
    state.constraint_residual = (state.entropy_j_kg_k - target).abs() / target.abs().max(1.0);
    Ok(state)
}

impl Solver {
    /// Solve chamber/throat/supersonic exit; pressures Pa, velocities m/s.
    pub fn rocket(
        &self,
        weights: &[f64],
        temperatures: &[f64],
        pc: f64,
        area_ratio: f64,
        flow: FlowMode,
    ) -> Result<RocketResult> {
        if !area_ratio.is_finite() || area_ratio < 1.0 || !pc.is_finite() || pc <= 0.0 {
            return Err(Error::Input(
                "Rocket requires finite positive pressure and area ratio >= 1".into(),
            ));
        }
        let h0 = self.reactant_enthalpy(weights, temperatures)?;
        let chamber_eq = self.solve(Mode::HP(h0, pc), weights)?;
        let chamber = if flow == FlowMode::FrozenAtChamber {
            frozen_state(&chamber_eq, chamber_eq.temperature_k, pc)?
        } else {
            chamber_eq
        };
        // Keep chamber condensates as initial guesses; unstable phases can still leave.
        let mut nozzle_solver = self.clone();
        nozzle_solver.insert.extend(
            chamber
                .amounts
                .iter()
                .filter(|(n, sp)| *n > 0.0 && sp.data().phase() != 0)
                .map(|(_, sp)| *sp),
        );
        let expand = |p: f64| {
            if flow == FlowMode::FrozenAtChamber {
                frozen_expansion(&chamber, p)
            } else {
                nozzle_solver.solve(Mode::SP(chamber.entropy_j_kg_k, p), weights)
            }
        };
        let sonic_residual = |ln_p: f64| {
            let st = expand(ln_p.exp())?;
            let speed2 = 2.0 * (chamber.enthalpy_j_kg - st.enthalpy_j_kg);
            Ok(speed2 / st.sound_speed().powi(2) - 1.0)
        };
        let mut low = pc * 0.5;
        for _ in 0..20 {
            if sonic_residual(low.ln())? > 0.0 {
                break;
            }
            low *= 0.5;
        }
        let pt = bracketed_root(low.ln(), (pc * 0.99999).ln(), sonic_residual)?.exp();
        let throat = expand(pt)?;
        let vt = (2.0 * (chamber.enthalpy_j_kg - throat.enthalpy_j_kg)).sqrt();
        let flux = throat.density_kg_m3 * vt;
        if !flux.is_finite() || flux <= 0.0 {
            return Err(Error::Convergence("Invalid throat mass flux".into()));
        }
        let exit_at = |p: f64| {
            if flow == FlowMode::FrozenAtThroat {
                frozen_expansion(&throat, p)
            } else {
                expand(p)
            }
        };
        let area_residual = |ln_p: f64| {
            let st = exit_at(ln_p.exp())?;
            let v = (2.0 * (chamber.enthalpy_j_kg - st.enthalpy_j_kg)).sqrt();
            Ok((flux / (st.density_kg_m3 * v) / area_ratio).ln())
        };
        let pe = if area_ratio == 1.0 {
            pt
        } else {
            let mut low = pt * 0.5;
            for _ in 0..40 {
                if area_residual(low.ln())? > 0.0 {
                    break;
                }
                low *= 0.5;
            }
            bracketed_root(low.ln(), pt.ln(), area_residual)?.exp()
        };
        let exit = exit_at(pe)?;
        let ve = (2.0 * (chamber.enthalpy_j_kg - exit.enthalpy_j_kg)).sqrt();
        let cstar = pc / flux;
        let vacuum_velocity = ve + pe * area_ratio / flux;
        if !vacuum_velocity.is_finite()
            || (flux / (exit.density_kg_m3 * ve) / area_ratio - 1.0).abs() > 2e-6
            || (vt / throat.sound_speed() - 1.0).abs() > 2e-6
        {
            return Err(Error::Convergence(
                "Nozzle failed sonic, mass-flow or energy closure".into(),
            ));
        }
        Ok(RocketResult {
            chamber: Station {
                state: chamber,
                velocity_m_s: 0.0,
                mach: 0.0,
                area_ratio: None,
            },
            throat: Station {
                mach: vt / throat.sound_speed(),
                state: throat,
                velocity_m_s: vt,
                area_ratio: Some(1.0),
            },
            exit: Station {
                mach: ve / exit.sound_speed(),
                state: exit,
                velocity_m_s: ve,
                area_ratio: Some(area_ratio),
            },
            cstar_m_s: cstar,
            vacuum_velocity_m_s: vacuum_velocity,
        })
    }
}

#[cfg(test)]
mod root_tests {
    use super::*;

    #[test]
    fn smooth_pressure_root_uses_few_calls_and_stays_bracketed() {
        let mut calls = 0;
        let result = bracketed_root(0.0, 10.0, |x| {
            assert!((0.0..=10.0).contains(&x));
            calls += 1;
            Ok(x.exp() - 20.0)
        })
        .unwrap();
        assert!((result - 20.0_f64.ln()).abs() < 1e-9);
        assert!(calls < 25, "{calls} expensive residual evaluations");
    }

    #[test]
    fn endpoint_roots_and_invalid_brackets_are_explicit() {
        assert_eq!(bracketed_root(0.0, 1.0, Ok).unwrap(), 0.0);
        assert_eq!(bracketed_root(0.0, 1.0, |x| Ok(x - 1.0)).unwrap(), 1.0);
        assert!(bracketed_root(0.0, 1.0, |x| Ok(x + 1.0)).is_err());
        assert!(bracketed_root(0.0, 1.0, |_| Ok(f64::NAN)).is_err());
        assert!(bracketed_root(0.0, 1.0, |_| Ok(f64::INFINITY)).is_err());
    }

    #[test]
    fn strongly_skewed_root_does_not_stagnate() {
        let root = bracketed_root(0.0, 2.0, |x| Ok(x.powi(10) - 1.0)).unwrap();
        assert!((root - 1.0).abs() < 1e-9);
    }
}
