//! Rocket chamber and nozzle calculations derived from NASA CEA `rocket.f90`.
//!
//! Apache-2.0 derived work, NASA CEA commit 2f79a647. Public quantities use SI.
use crate::{
    Error, R, Result,
    solver::{Mode, Solver, State},
    transport::Transport,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowMode {
    Equilibrium,
    FrozenAtChamber,
    FrozenAtThroat,
}

/// Exit schedules run in pressure, subsonic-area, supersonic-area order.
/// `n_frz` uses NASA's one-based station indices; zero means equilibrium.
#[derive(Debug, Clone)]
pub struct RocketOptions<'a> {
    pub pi_p: &'a [f64],
    pub subar: &'a [f64],
    pub supar: &'a [f64],
    pub iac: bool,
    pub ac_at: Option<f64>,
    /// Mass flow per chamber area, kg/(s m²).
    pub mdot: Option<f64>,
    pub n_frz: usize,
    /// Assigned stagnation enthalpy, J/kg; overrides reactant-stream enthalpy.
    pub hc: Option<f64>,
    /// Assigned stagnation temperature, K; mutually exclusive with hc.
    pub tc: Option<f64>,
    /// Initial temperature for the assigned-enthalpy chamber iteration, K.
    pub tc_est: Option<f64>,
    /// Fill `Station::transport` at every station, as NASA's `transport` flag.
    pub transport: bool,
}
impl Default for RocketOptions<'_> {
    fn default() -> Self {
        Self {
            pi_p: &[],
            subar: &[],
            supar: &[],
            iac: true,
            ac_at: None,
            mdot: None,
            n_frz: 0,
            hc: None,
            tc: None,
            tc_est: None,
            transport: false,
        }
    }
}

#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct Station {
    pub name: &'static str,
    pub state: State,
    pub velocity_m_s: f64,
    pub sound_speed_m_s: f64,
    pub mach: f64,
    pub area_ratio: Option<f64>,
    pub cstar_m_s: f64,
    pub thrust_coefficient: f64,
    /// Ambient-matched specific impulse expressed as velocity, m/s.
    pub isp_m_s: f64,
    pub isp_vacuum_m_s: f64,
    /// `state.transport()`, set only when `RocketOptions::transport` is true.
    /// As in NASA, a frozen station uses its frozen composition and gas Cp.
    pub transport: Option<Transport>,
}
#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct RocketResult {
    /// IAC infinity station or FAC combustor-end station.
    pub chamber: Station,
    pub throat: Station,
    /// Last completed exit, or throat when the schedule contains no exits.
    pub exit: Station,
    pub cstar_m_s: f64,
    pub vacuum_velocity_m_s: f64,
    /// FAC includes injector, infinity, combustor end and throat before exits.
    pub stations: Vec<Station>,
    pub throat_index: usize,
    /// False when NASA's frozen fit limits stop an expansion early.
    pub complete: bool,
    pub stop_reason: Option<String>,
}

fn frozen_bounds(reference: &State) -> (f64, f64) {
    let mut gas_min = f64::INFINITY;
    let mut low = 0.0_f64;
    let mut high = f64::INFINITY;
    for &(n, sp) in &reference.amounts {
        if n <= 1e-12 {
            continue;
        }
        let d = sp.data();
        if d.phase() == 0 {
            gas_min = gas_min.min(d.valid_temperature_range().0);
        } else {
            // NASA freezes a condensed phase in its reference fit interval.
            let ranges = d.temperature_data();
            let interval = ranges
                .iter()
                .rev()
                .find(|r| reference.temperature_k > r.0)
                .unwrap_or(&ranges[0]);
            low = low.max(interval.0 - 50.0);
            high = high.min(interval.1 + 50.0);
        }
    }
    (low.max(0.8 * gas_min), high)
}

fn frozen_expansion(reference: &State, p: f64) -> Result<State> {
    let (min_t, max_t) = frozen_bounds(reference);
    if reference.frozen(min_t, p)?.entropy_j_kg_k > reference.entropy_j_kg_k
        || (max_t.is_finite()
            && reference.frozen(max_t, p)?.entropy_j_kg_k < reference.entropy_j_kg_k)
    {
        return Err(Error::Domain(format!(
            "Frozen calculations stopped: entropy requires a temperature outside {min_t}..{max_t} K"
        )));
    }
    let mut t = reference.temperature_k
        * (p / reference.pressure_pa).powf((reference.gamma_s - 1.0) / reference.gamma_s);
    t = t.clamp(min_t, max_t);
    for _ in 0..30 {
        let mut st = reference.frozen(t, p)?;
        let ds = reference.entropy_j_kg_k - st.entropy_j_kg_k;
        if ds.abs() / reference.entropy_j_kg_k.abs().max(1.0) < 1e-11 {
            if t < min_t || t > max_t {
                return Err(Error::Domain(format!(
                    "Frozen calculations stopped at {t:.6} K: allowed reference fit range {min_t}..{max_t} K"
                )));
            }
            st.constraint_residual = ds.abs() / reference.entropy_j_kg_k.abs().max(1.0);
            return Ok(st);
        }
        t = (t * (ds / st.cp_frozen_j_kg_k).clamp(-0.5, 0.5).exp()).clamp(min_t, max_t);
    }
    Err(Error::Convergence(
        "Frozen temperature did not converge".into(),
    ))
}

fn station(name: &'static str, state: State, h0: f64, flux: Option<f64>) -> Result<Station> {
    let speed2 = if flux.is_some() {
        2.0 * (h0 - state.enthalpy_j_kg)
    } else {
        0.0
    };
    if !speed2.is_finite() || speed2 < 0.0 {
        return Err(Error::Convergence(
            "Nozzle has negative kinetic energy".into(),
        ));
    }
    let velocity = speed2.sqrt();
    let area = flux.map(|g| g / (state.density_kg_m3 * velocity));
    if area.is_some_and(|a| !a.is_finite() || a <= 0.0) {
        return Err(Error::Convergence("Invalid nozzle mass flux".into()));
    }
    Ok(Station {
        name,
        sound_speed_m_s: state.sound_speed(),
        mach: velocity / state.sound_speed(),
        state,
        velocity_m_s: velocity,
        area_ratio: area,
        cstar_m_s: 0.0,
        thrust_coefficient: 0.0,
        isp_m_s: velocity,
        isp_vacuum_m_s: 0.0,
        transport: None,
    })
}

/// An isentrope continues from the previous equilibrium solution.
struct Nozzle<'a> {
    solver: &'a Solver,
    start: Option<crate::equilibrium_core::EquilibriumSeed>,
    weights: &'a [f64],
    reference: &'a State,
    frozen: bool,
    h0: f64,
}
impl Nozzle<'_> {
    fn state(&mut self, p: f64) -> Result<State> {
        if self.frozen {
            return frozen_expansion(self.reference, p);
        }
        let st = self.solver.solve_with_start(
            Mode::SP(self.reference.entropy_j_kg_k, p),
            self.weights,
            self.start.as_ref(),
        )?;
        self.start.clone_from(&st.equilibrium_seed);
        Ok(st)
    }

    fn throat(&mut self, p0: f64) -> Result<Station> {
        let gamma = self.reference.gamma_s;
        let mut p = p0 / ((gamma + 1.0) / 2.0).powf(gamma / (gamma - 1.0));
        for _ in 0..40 {
            let st = self.state(p)?;
            let speed2 = 2.0 * (self.h0 - st.enthalpy_j_kg);
            let mach2 = speed2 / st.sound_speed().powi(2);
            if !mach2.is_finite() || mach2 <= 0.0 {
                return Err(Error::Convergence("Nonphysical throat iterate".into()));
            }
            if (mach2 - 1.0).abs() < 2e-8 {
                let flux = st.density_kg_m3 * speed2.sqrt();
                return station("throat", st, self.h0, Some(flux));
            }
            // NASA RP-1311 Eq. 6.17, using the current isentropic exponent.
            p *= (1.0 + st.gamma_s * mach2) / (1.0 + st.gamma_s);
            p = p.clamp(p0 * 1e-6, p0 * (1.0 - 1e-8));
        }
        Err(Error::Convergence("Sonic throat did not converge".into()))
    }

    fn area(&mut self, p0: f64, throat: &Station, area: f64, supersonic: bool) -> Result<Station> {
        if area == 1.0 {
            let mut st = throat.clone();
            st.name = "exit";
            return Ok(st);
        }
        let flux = throat.state.density_kg_m3 * throat.velocity_m_s;
        let xt = (p0 / throat.state.pressure_pa).ln();
        let loga = area.ln();
        let mut x = if supersonic {
            if area < 2.0 {
                xt + (3.294 * loga.powi(2) + 1.535 * loga).sqrt()
            } else {
                throat.state.gamma_s + 1.4 * loga
            }
        } else {
            let mut x = xt / (area + 10.587 * loga.powi(3) + 9.454 * loga);
            if area < 1.09 {
                x *= 0.9;
            } else if area > 10.0 {
                x /= area;
            }
            x
        };
        let mut lo = if supersonic { xt } else { 0.0 };
        let mut hi = if supersonic { f64::INFINITY } else { xt };
        for _ in 0..40 {
            let st = station("exit", self.state(p0 * (-x).exp())?, self.h0, Some(flux))?;
            let actual = st
                .area_ratio
                .ok_or_else(|| Error::Convergence("Missing nozzle area".into()))?;
            let residual = (area / actual).ln();
            let m2 = st.mach.powi(2);
            let update = st.state.gamma_s * m2 / (m2 - 1.0) * residual;
            if residual.abs() < 2e-8
                || (!supersonic
                    && (actual / area - 1.0).abs() < 5e-4
                    && ((actual / area - 1.0).abs() <= 4e-5 || update.abs() < 4e-5))
            {
                return Ok(st);
            }
            if (actual < area) == supersonic {
                lo = x;
            } else {
                hi = x;
            }
            // NASA RP-1311 Eq. 6.23; safeguard the two area branches near Mach 1.
            let next = x + update;
            x = if next.is_finite() && next > lo && next < hi {
                next
            } else if hi.is_finite() {
                (lo + hi) / 2.0
            } else {
                x + 1.0
            };
        }
        Err(Error::Convergence(
            "Assigned nozzle area did not converge".into(),
        ))
    }
}

impl Solver {
    /// Backward-compatible IAC chamber/throat/supersonic exit solve.
    pub fn rocket(
        &self,
        weights: &[f64],
        temperatures: &[f64],
        pc: f64,
        area_ratio: f64,
        flow: FlowMode,
    ) -> Result<RocketResult> {
        let result = self.rocket_stations(
            weights,
            temperatures,
            pc,
            &RocketOptions {
                supar: &[area_ratio],
                n_frz: match flow {
                    FlowMode::Equilibrium => 0,
                    FlowMode::FrozenAtChamber => 1,
                    FlowMode::FrozenAtThroat => 2,
                },
                ..RocketOptions::default()
            },
        )?;
        if !result.complete {
            return Err(Error::Domain(
                result
                    .stop_reason
                    .unwrap_or_else(|| "Frozen exit was not completed".into()),
            ));
        }
        Ok(result)
    }

    /// Solve IAC or FAC and any number of pressure/area stations, in SI units.
    pub fn rocket_stations(
        &self,
        weights: &[f64],
        temperatures: &[f64],
        pc: f64,
        options: &RocketOptions<'_>,
    ) -> Result<RocketResult> {
        let o = options;
        self.feeds(weights)?;
        if !pc.is_finite()
            || pc <= 0.0
            || o.pi_p.iter().any(|p| !p.is_finite() || *p <= 1.0)
            || o.subar.iter().any(|a| !a.is_finite() || *a <= 1.0)
            || o.supar.iter().any(|a| !a.is_finite() || *a < 1.0)
            || o.hc.is_some_and(|h| !h.is_finite())
            || o.tc.is_some_and(|t| !t.is_finite() || t <= 0.0)
            || o.tc_est.is_some_and(|t| !t.is_finite() || t <= 0.0)
            || (o.hc.is_some() && o.tc.is_some())
        {
            return Err(Error::Input(
                "Invalid rocket pressure, area or assigned chamber state".into(),
            ));
        }
        if (o.iac && (o.ac_at.is_some() || o.mdot.is_some()))
            || (!o.iac && (o.ac_at.is_some() == o.mdot.is_some()))
            || o.ac_at.is_some_and(|a| !a.is_finite() || a < 1.0)
            || o.mdot.is_some_and(|m| !m.is_finite() || m <= 0.0)
        {
            return Err(Error::Input(
                "FAC requires exactly one of ac_at >= 1 or mdot > 0; IAC takes neither".into(),
            ));
        }
        let throat_index = if o.iac { 1 } else { 3 };
        let num_pts = throat_index + 1 + o.pi_p.len() + o.subar.len() + o.supar.len();
        if o.n_frz > num_pts {
            return Err(Error::Input(
                "Freeze index exceeds the requested station count".into(),
            ));
        }
        let assigned_h = match (o.hc, o.tc) {
            (Some(h), _) => h,
            (_, Some(_)) => 0.0,
            _ => self.reactant_enthalpy(weights, temperatures)?,
        };
        let chamber_at = |p| {
            if let Some(t) = o.tc {
                return self.solve(Mode::TP(t, p), weights);
            }
            if let Some(mut t) = o.tc_est {
                // Public solve has no temperature seed; honor tc_est through TP heat-capacity updates.
                for _ in 0..40 {
                    let st = self.solve(Mode::TP(t, p), weights)?;
                    let dh = assigned_h - st.enthalpy_j_kg;
                    let scale = assigned_h
                        .abs()
                        .max(R * t / st.molecular_weight_kg_kmol)
                        .max(1.0);
                    if dh.abs() < 1e-9 * scale {
                        return Ok(st);
                    }
                    t += (dh / st.cp_equilibrium_j_kg_k).clamp(-0.25 * t, 0.25 * t);
                }
                return Err(Error::Convergence(
                    "Estimated chamber enthalpy iteration did not converge".into(),
                ));
            }
            self.solve(Mode::HP(assigned_h, p), weights)
        };
        let mut injector = chamber_at(pc)?;
        let h0 = injector.enthalpy_j_kg;
        let mut stations;
        let throat;
        let infinity;
        if o.iac {
            if o.n_frz == 1 {
                injector = injector.frozen(injector.temperature_k, pc)?;
            }
            infinity = injector.clone();
            let mut nozzle = Nozzle {
                solver: self,
                weights,
                start: infinity.equilibrium_seed.clone(),
                reference: &infinity,
                frozen: o.n_frz == 1,
                h0,
            };
            throat = nozzle.throat(pc)?;
            stations = vec![station("infinity", injector, h0, None)?, throat.clone()];
        } else {
            let mut ac = o.ac_at.unwrap_or(2.0);
            let ratio = |a: f64| (1.0257 - 1.2318 * a) / (1.0 - 1.26505 * a);
            let mut pinf = pc * ratio(ac);
            let mut closed = None;
            for _ in 0..100 {
                let inf = if o.n_frz == 1 {
                    injector.frozen(injector.temperature_k, pinf)?
                } else {
                    chamber_at(pinf)?
                };
                let mut nozzle = Nozzle {
                    solver: self,
                    weights,
                    start: inf.equilibrium_seed.clone(),
                    reference: &inf,
                    frozen: (1..=2).contains(&o.n_frz),
                    h0,
                };
                let mut thr = nozzle.throat(pinf)?;
                let mut consistent_combustor = None;
                for _ in 0..40 {
                    let flux = thr.state.density_kg_m3 * thr.velocity_m_s;
                    if let Some(mdot) = o.mdot {
                        ac = flux / mdot;
                    }
                    if ac < 1.0 || !ac.is_finite() {
                        return Err(Error::Input(
                            "FAC mass flux exceeds the sonic chamber limit".into(),
                        ));
                    }
                    let mut combustor = nozzle.area(pinf, &thr, ac, false)?;
                    combustor.name = "combustor";
                    if o.n_frz == 3 {
                        let mut frozen = Nozzle {
                            solver: self,
                            weights,
                            start: combustor.state.equilibrium_seed.clone(),
                            reference: &combustor.state,
                            frozen: true,
                            h0,
                        };
                        let frozen_throat = frozen.throat(pinf)?;
                        let frozen_flux =
                            frozen_throat.state.density_kg_m3 * frozen_throat.velocity_m_s;
                        thr = frozen_throat;
                        if (frozen_flux / flux - 1.0).abs() > 2e-8 {
                            continue;
                        }
                    }
                    consistent_combustor = Some(combustor);
                    break;
                }
                let combustor = consistent_combustor.ok_or_else(|| {
                    Error::Convergence(
                        "Frozen FAC combustor/throat mass closure did not converge".into(),
                    )
                })?;
                let check = combustor.state.pressure_pa
                    + combustor.state.density_kg_m3 * combustor.velocity_m_s.powi(2);
                if (check / pc - 1.0).abs() < 2e-8 {
                    closed = Some((inf, combustor, thr));
                    break;
                }
                pinf *= pc / check;
            }
            let (inf, combustor, thr) = closed.ok_or_else(|| {
                Error::Convergence("FAC injector momentum closure did not converge".into())
            })?;
            infinity = inf;
            throat = thr;
            stations = vec![
                station("injector", injector, h0, None)?,
                Station {
                    // NASA leaves this auxiliary FAC station's sonic output at zero.
                    sound_speed_m_s: 0.0,
                    ..station("infinity", infinity.clone(), h0, None)?
                },
                combustor,
                throat.clone(),
            ];
        }
        let flux = throat.state.density_kg_m3 * throat.velocity_m_s;
        if !flux.is_finite() || flux <= 0.0 || (throat.mach - 1.0).abs() > 2e-6 {
            return Err(Error::Convergence(
                "Nozzle failed sonic or mass-flow closure".into(),
            ));
        }
        let mut stop_reason = None;
        for (kind, schedule) in [(0, o.pi_p), (1, o.subar), (2, o.supar)] {
            // NASA omits subsonic schedules when freezing after station 1.
            if kind == 1 && o.n_frz > 1 {
                continue;
            }
            for &value in schedule {
                let frozen = o.n_frz > 0 && stations.len() >= o.n_frz;
                let reference = if frozen && !o.iac && o.n_frz == 1 {
                    &infinity
                } else if frozen {
                    &stations[o.n_frz - 1].state
                } else {
                    &infinity
                };
                if frozen && kind == 0 && value < pc / reference.pressure_pa {
                    continue;
                }
                if frozen
                    && kind == 2
                    && o.n_frz >= 3
                    && value <= stations[o.n_frz - 1].area_ratio.unwrap_or(0.0)
                {
                    continue;
                }
                let mut nozzle = Nozzle {
                    solver: self,
                    weights,
                    start: stations
                        .last()
                        .and_then(|st| st.state.equilibrium_seed.clone()),
                    reference,
                    frozen,
                    h0,
                };
                let result = if kind == 0 {
                    nozzle
                        .state(pc / value)
                        .and_then(|st| station("exit", st, h0, Some(flux)))
                } else {
                    // Native n_frz=1 keeps subar in equilibrium, with the frozen throat flux.
                    if kind == 1 {
                        nozzle.frozen = false;
                        nozzle.reference = &infinity;
                    }
                    nozzle.area(infinity.pressure_pa, &throat, value, kind == 2)
                };
                let st = match result {
                    Ok(st) => st,
                    Err(Error::Domain(message)) if nozzle.frozen => {
                        stop_reason = Some(message);
                        break;
                    }
                    Err(err) => return Err(err),
                };
                let area_tolerance = if kind == 1 { 5e-4 } else { 2e-6 };
                if kind != 0 && (st.area_ratio.unwrap_or(0.0) / value - 1.0).abs() > area_tolerance
                {
                    return Err(Error::Convergence(
                        "Nozzle failed assigned mass-flow closure".into(),
                    ));
                }
                stations.push(st);
            }
            if stop_reason.is_some() {
                break;
            }
        }
        let cstar = infinity.pressure_pa / flux;
        for st in &mut stations {
            st.cstar_m_s = cstar;
            if let Some(area) = st.area_ratio {
                st.thrust_coefficient = st.velocity_m_s / cstar;
                if (st.state.density_kg_m3 * st.velocity_m_s * area / flux - 1.0).abs() > 2e-6 {
                    return Err(Error::Convergence(
                        "Nozzle failed station mass-flow closure".into(),
                    ));
                }
                st.isp_vacuum_m_s = st.velocity_m_s
                    + st.state.pressure_pa / (st.state.density_kg_m3 * st.velocity_m_s);
                let energy = st.state.enthalpy_j_kg + st.velocity_m_s.powi(2) / 2.0;
                if !st.isp_vacuum_m_s.is_finite() || (energy - h0).abs() > 2e-6 * h0.abs().max(1.0)
                {
                    return Err(Error::Convergence("Nozzle failed energy closure".into()));
                }
            }
            if o.transport {
                st.transport = Some(st.state.transport()?);
            }
        }
        let exit = stations
            .last()
            .ok_or_else(|| Error::Convergence("Empty rocket result".into()))?
            .clone();
        Ok(RocketResult {
            chamber: stations[if o.iac { 0 } else { 2 }].clone(),
            throat: stations[throat_index].clone(),
            vacuum_velocity_m_s: exit.isp_vacuum_m_s,
            exit,
            cstar_m_s: cstar,
            stations,
            throat_index,
            complete: stop_reason.is_none(),
            stop_reason,
        })
    }
}
