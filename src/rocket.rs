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

fn frozen_expansion(reference: &State, p: f64, mut t: f64) -> Result<State> {
    let (min_t, max_t) = frozen_bounds(reference);
    let mut converged = false;
    for _ in 0..8 {
        let mut st = reference.frozen(t, p)?;
        if converged {
            if t < min_t || t > max_t {
                return Err(Error::Domain(format!(
                    "Frozen calculations stopped: temperature outside {min_t}..{max_t} K"
                )));
            }
            st.entropy_j_kg_k = reference.entropy_j_kg_k;
            st.gibbs_energy_j_kg = st.enthalpy_j_kg - t * st.entropy_j_kg_k;
            st.cv_frozen_j_kg_k = 0.0;
            st.cv_equilibrium_j_kg_k = 0.0;
            return Ok(st);
        }
        let dlnt = (reference.entropy_j_kg_k - st.entropy_j_kg_k) / st.cp_frozen_j_kg_k;
        t = (t.ln() + dlnt).exp();
        if dlnt.abs() < 0.5e-4 {
            converged = true;
        }
    }
    Err(Error::Domain(
        "Frozen calculations did not converge in 8 iterations".into(),
    ))
}

fn station(name: &'static str, state: State, h0: f64, flux: Option<f64>) -> Result<Station> {
    let speed2 = if flux.is_some() {
        2.0 * (h0 - state.iteration_enthalpy())
    } else {
        0.0
    };
    let velocity = speed2.sqrt();
    let area = flux
        .map(|g| g * state.gas_moles * R * state.temperature_k / (state.pressure_pa * velocity));
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
    temperature: f64,
    h0: f64,
}
impl Nozzle<'_> {
    fn state(&mut self, p: f64) -> Result<State> {
        if self.frozen {
            let st = frozen_expansion(self.reference, p, self.temperature)?;
            self.temperature = st.temperature_k;
            return Ok(st);
        }
        let mut st = self.solver.solve_with_start(
            Mode::SP(self.reference.entropy_j_kg_k, p),
            self.weights,
            self.start.as_ref(),
        )?;
        if st.converged {
            self.start.clone_from(&st.equilibrium_seed);
        }
        st.pressure_pa = p;
        Ok(st)
    }

    fn throat(&mut self, p0: f64, gamma: f64) -> Result<Station> {
        let mut p = p0 / ((gamma + 1.0) / 2.0).powf(gamma / (gamma - 1.0));
        let mut last = None;
        let mut melt_t = 0.0_f64;
        for i in 1..=22 {
            let st = self.state(p)?;
            let speed2 = 2.0 * (self.h0 - st.iteration_enthalpy());
            let mach2 = speed2 / st.sound_speed().powi(2);
            let mut current = station("throat", st, self.h0, Some(1.0))?;
            let flux = current.state.pressure_pa * current.velocity_m_s
                / (current.state.gas_moles * R * current.state.temperature_k);
            current.area_ratio = Some(1.0);
            if self.frozen {
                current.state.cv_frozen_j_kg_k = self.reference.cv_frozen_j_kg_k;
                current.state.cv_equilibrium_j_kg_k = self.reference.cv_equilibrium_j_kg_k;
            }
            if (mach2 - 1.0).abs() / mach2 <= 0.4e-4 {
                return Ok(current);
            }
            let next = p * (1.0 + current.state.gamma_s * mach2) / (1.0 + current.state.gamma_s);
            p = if i <= 3 {
                if !self.frozen
                    && current
                        .state
                        .equilibrium_seed
                        .as_ref()
                        .is_some_and(|s| s.j_sol.is_some())
                {
                    melt_t = current.state.temperature_k;
                    next
                } else if !self.frozen && melt_t.abs() >= 1e-6 {
                    current.state.pressure_pa = p
                        * ((melt_t / current.state.temperature_k).ln()
                            * current.state.cp_frozen_j_kg_k
                            / (current.state.gas_moles * R * current.state.dlnv_dlnt))
                            .exp();
                    return Ok(current);
                } else {
                    next
                }
            } else {
                next.max(p) - (next - p).abs() / 20.0
            };
            last = Some((current, flux));
        }
        let mut current = last.unwrap().0;
        current.state.pressure_pa = p;
        Ok(current)
    }

    fn area(
        &mut self,
        p0: f64,
        throat: &Station,
        area: f64,
        supersonic: bool,
        gamma_estimate: f64,
    ) -> Result<Station> {
        let flux = throat.state.pressure_pa * throat.velocity_m_s
            / (throat.state.gas_moles * R * throat.state.temperature_k);
        let xt = (p0 / throat.state.pressure_pa).ln();
        let loga = area.ln();
        let mut x = if supersonic {
            if area < 2.0 {
                xt + (3.294 * loga.powi(2) + 1.535 * loga).sqrt()
            } else {
                gamma_estimate + 1.4 * loga
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
        let mut last = None;
        for _ in 0..10 {
            let st = station("exit", self.state(p0 / x.exp())?, self.h0, Some(flux))?;
            let actual = st.area_ratio.unwrap();
            let speed2 = st.velocity_m_s.powi(2);
            let sound2 = st.sound_speed_m_s.powi(2);
            let update = st.state.gamma_s * speed2 / (speed2 - sound2) * (area.ln() - actual.ln());
            if (actual - area).abs() / area <= 4e-5 || update.abs() < 4e-5 {
                return Ok(st);
            }
            x += update;
            last = Some(st);
        }
        Ok(last.unwrap())
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
            || o.supar.iter().any(|a| !a.is_finite() || *a <= 1.0)
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
            if let Some(t) = o.tc_est {
                let seed = self.initial_seed(weights, t, None)?;
                return self.solve_with_start(Mode::HP(assigned_h, p), weights, Some(&seed));
            }
            self.solve(Mode::HP(assigned_h, p), weights)
        };
        let mut injector = chamber_at(pc)?;
        let h0 = if o.tc.is_some() {
            injector.iteration_enthalpy()
        } else {
            assigned_h
        };
        let mut stations;
        let throat;
        let infinity;
        if o.iac {
            if o.n_frz == 1 {
                injector.gamma_s = injector.cp_frozen_j_kg_k
                    / (injector.cp_frozen_j_kg_k - injector.gas_moles * R);
                injector.equilibrium = false;
            }
            infinity = injector.clone();
            let mut nozzle = Nozzle {
                solver: self,
                weights,
                start: infinity.equilibrium_seed.clone(),
                reference: &infinity,
                frozen: o.n_frz == 1,
                temperature: infinity.temperature_k,
                h0,
            };
            throat = nozzle.throat(pc, infinity.gamma_s)?;
            stations = vec![station("infinity", injector, h0, None)?, throat.clone()];
        } else {
            let mut ac = o.ac_at.unwrap_or(2.0);
            let ratio = |a: f64| (1.0257 - 1.2318 * a) / (1.0 - 1.26505 * a);
            let mut prat = 0.0;
            let mut pinf = pc;
            for _ in 0..4 {
                prat = ratio(ac);
                pinf = pc * prat;
                if o.ac_at.is_some() {
                    break;
                }
                ac = pinf / (o.mdot.unwrap() * 2350.0);
                if ac < 1.0 {
                    return Err(Error::Input(
                        "FAC mass flow gives a contraction ratio below one".into(),
                    ));
                }
            }
            let mut inf_seed = self.initial_seed(weights, o.tc_est.unwrap_or(3800.0), None)?;
            let mut closed = None;
            for _ in 0..100 {
                let mode = if let Some(t) = o.tc {
                    Mode::TP(t, pinf)
                } else {
                    Mode::HP(assigned_h, pinf)
                };
                let inf = self.solve_with_start(mode, weights, Some(&inf_seed))?;
                if inf.converged {
                    inf_seed = inf.equilibrium_seed.clone().unwrap();
                }
                let mut nozzle = Nozzle {
                    solver: self,
                    weights,
                    start: inf.equilibrium_seed.clone(),
                    reference: &inf,
                    frozen: false,
                    temperature: inf.temperature_k,
                    h0,
                };
                let thr = nozzle.throat(pinf, injector.gamma_s)?;
                nozzle.start = inf.equilibrium_seed.clone();
                let mut combustor = nozzle.area(pinf, &thr, ac, false, inf.gamma_s)?;
                combustor.name = "combustor";
                let check = combustor.state.pressure_pa
                    + combustor.velocity_m_s.powi(2) * combustor.state.pressure_pa
                        / (R * combustor.state.gas_moles * combustor.state.temperature_k);
                let done = if o.ac_at.is_some() {
                    let done = (pc - check).abs() / pc <= 2e-5;
                    if !done {
                        pinf *= pc / check;
                    }
                    done
                } else {
                    let flux = thr.state.pressure_pa * thr.velocity_m_s
                        / (R * thr.state.gas_moles * thr.state.temperature_k);
                    ac = flux / o.mdot.unwrap();
                    let mut previous = prat;
                    prat = ratio(ac);
                    let done = (pc - check).abs() / check <= 2e-5;
                    if !done {
                        let mut pjrat = check / pc;
                        for _ in 0..2 {
                            let pr = pjrat * previous / prat;
                            pinf /= pr;
                            ac /= pr;
                            combustor.area_ratio = Some(ac);
                            previous = prat;
                            pjrat = 1.0;
                            prat = ratio(ac);
                        }
                    }
                    done
                };
                closed = Some((inf, combustor, thr));
                if done {
                    break;
                }
            }
            let (inf, combustor, mut thr) = closed.unwrap();
            // NASA leaves the FAC infinity/combustor gamma workspace at zero here.
            if (1..=3).contains(&o.n_frz) {
                let reference = match o.n_frz {
                    1 => &injector,
                    2 => &inf,
                    _ => &combustor.state,
                };
                let mut nozzle = Nozzle {
                    solver: self,
                    weights,
                    start: None,
                    reference,
                    frozen: true,
                    temperature: reference.temperature_k,
                    h0,
                };
                thr = nozzle.throat(pinf, if o.n_frz == 1 { reference.gamma_s } else { 0.0 })?;
            }
            infinity = inf;
            throat = thr;
            stations = vec![
                station("injector", injector, h0, None)?,
                Station {
                    sound_speed_m_s: 0.0,
                    ..station("infinity", infinity.clone(), h0, None)?
                },
                combustor,
                throat.clone(),
            ];
        }
        let flux = throat.state.pressure_pa * throat.sound_speed_m_s
            / (throat.state.gas_moles * R * throat.state.temperature_k);
        let mut stop_reason = None;
        for (kind, schedule) in [(0, o.pi_p), (1, o.subar), (2, o.supar)] {
            let mut frozen_temperature = None;
            // NASA omits subsonic schedules when freezing after station 1.
            if kind == 1 && o.n_frz > 1 {
                continue;
            }
            for &value in schedule {
                let frozen = o.n_frz > 0 && stations.len() >= o.n_frz;
                let reference = if frozen {
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
                    temperature: *frozen_temperature.get_or_insert_with(|| {
                        if kind == 0 {
                            infinity.temperature_k
                        } else {
                            stations.last().unwrap().state.temperature_k
                        }
                    }),
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
                    nozzle.area(
                        infinity.pressure_pa,
                        &throat,
                        value,
                        kind == 2,
                        if o.iac {
                            throat.state.gamma_s
                        } else {
                            infinity.gamma_s
                        },
                    )
                };
                let st = match result {
                    Ok(st) => st,
                    Err(Error::Domain(message)) if nozzle.frozen => {
                        stop_reason = Some(message);
                        break;
                    }
                    Err(err) => return Err(err),
                };
                stations.push(st);
            }
            if stop_reason.is_some() {
                break;
            }
        }
        let performance_h0 = stations[0].state.iteration_enthalpy();
        let throat_isp = (2.0 * (performance_h0 - throat.state.iteration_enthalpy())).sqrt();
        let throat_flux = throat.state.pressure_pa * throat_isp
            / (throat.state.gas_moles * R * throat.state.temperature_k);
        let cstar = infinity.pressure_pa / throat_flux;
        for st in &mut stations {
            st.cstar_m_s = cstar;
            if st.area_ratio.is_some() {
                st.isp_m_s = (2.0 * (performance_h0 - st.state.iteration_enthalpy())).sqrt();
                st.thrust_coefficient = st.isp_m_s / cstar;
                st.isp_vacuum_m_s =
                    st.isp_m_s + st.state.gas_moles * R * st.state.temperature_k / st.isp_m_s;
            }
            if o.transport {
                let tr = st.state.transport()?;
                st.state.cp_frozen_j_kg_k = tr.cp_frozen_j_kg_k;
                if st.state.equilibrium {
                    st.state.cv_frozen_j_kg_k = tr.cp_frozen_j_kg_k - st.state.gas_moles * R;
                }
                st.transport = Some(tr);
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
