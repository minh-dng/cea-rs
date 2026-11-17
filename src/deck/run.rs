//! Input schedules mapped to the existing public solver APIs, without changing numerics.
use super::{AmountBasis, Deck, Problem, Ratio, Role, invalid, parse_deck, positive};
use crate::{
    Error, Result,
    detonation::DetonationResult,
    rocket::{RocketOptions, RocketResult},
    shock::{ShockResult, ShockSpeed},
    solver::{Mode, Solver, State},
    transport::Transport,
};

#[derive(Debug, Clone, Copy)]
pub enum CaseInput {
    Equilibrium(Mode),
    Rocket {
        pressure_pa: f64,
        freeze_station: usize,
    },
    Shock {
        temperature_k: f64,
        pressure_pa: f64,
        speed: ShockSpeed,
        incident_frozen: bool,
        reflected_frozen: bool,
    },
    Detonation {
        temperature_k: f64,
        pressure_pa: f64,
    },
}

#[derive(Debug, Clone)]
pub struct Case<'a> {
    pub deck: &'a Deck,
    /// One-based case number within this deck, in NASA report order.
    pub index: usize,
    /// Reactant mass weights in input order (normalization is solver-owned).
    pub weights: Vec<f64>,
    pub input: CaseInput,
}

#[derive(Debug, Clone)]
pub enum CaseOutput {
    Equilibrium(Box<State>),
    Rocket(Box<RocketResult>),
    Shock(Box<ShockResult>),
    Detonation(Box<DetonationResult>),
}

#[cfg_attr(feature = "python", pyo3::pyclass(frozen, skip_from_py_object))]
#[derive(Debug, Clone)]
pub struct DeckResult {
    pub case_name: String,
    pub case_index: usize,
    pub weights: Vec<f64>,
    pub input: CaseInput,
    pub output: CaseOutput,
    /// Equilibrium/burned state; all rocket stations; or shock stations 1, 2, [5].
    /// Empty unless `output transport` was requested.
    pub transport: Vec<Transport>,
    /// Harmonic mass-weighted inlet density, available when every density was supplied.
    pub reactant_density_kg_m3: Option<f64>,
}

impl Deck {
    fn solver(&self) -> Result<Solver> {
        Solver::new(
            self.reactants.iter().map(|r| r.reactant.clone()).collect(),
            true,
            self.ions,
            self.only.clone(),
            self.insert.clone(),
        )?
        .with_omit(self.omit.clone())?
        .with_trace(self.output.trace)
    }

    fn temperatures(&self) -> Vec<f64> {
        self.reactants
            .iter()
            .map(|r| r.temperature_k.unwrap_or(r.reactant.reference_temperature))
            .collect()
    }

    fn weights(&self, solver: &Solver) -> Result<Vec<Vec<f64>>> {
        let amounts: Vec<_> = self.reactants.iter().map(|r| r.amount).collect();
        let weights = if self.reactants[0].basis == AmountBasis::Moles {
            solver.moles_to_weights(&amounts)?
        } else {
            amounts
        };
        let Some((kind, ratios)) = &self.ratios else {
            return Ok(vec![weights]);
        };
        let ox: Vec<_> = self
            .reactants
            .iter()
            .zip(&weights)
            .map(|(r, w)| if r.role == Role::Oxidant { *w } else { 0.0 })
            .collect();
        let fu: Vec<_> = self
            .reactants
            .iter()
            .zip(&weights)
            .map(|(r, w)| if r.role == Role::Fuel { *w } else { 0.0 })
            .collect();
        ratios
            .iter()
            .map(|&ratio| {
                let of = match kind {
                    Ratio::OxidantFuel => ratio,
                    Ratio::FuelOxidant => 1.0 / positive(ratio, "F/O")?,
                    Ratio::PercentFuel => {
                        positive(ratio, "%fuel")?;
                        if ratio > 100.0 {
                            return Err(invalid("%fuel must not exceed 100"));
                        }
                        (100.0 - ratio) / ratio
                    }
                    Ratio::Phi => solver.weight_eq_ratio_to_of_ratio(&ox, &fu, ratio)?,
                    Ratio::Chemical => solver.chem_eq_ratio_to_of_ratio(&ox, &fu, ratio)?,
                };
                solver.of_ratio_to_weights(&ox, &fu, of)
            })
            .collect()
    }

    /// Expand schedules without running the application solvers. Individual cases
    /// can be run independently so a failed solve does not hide later diagnostics.
    pub fn cases(&self) -> Result<Vec<Case<'_>>> {
        if self.reactants.is_empty() {
            return Err(invalid("No reactants"));
        }
        if self.frozen && !matches!(self.problem, Problem::Rocket | Problem::Shock) {
            return Err(Error::Unsupported(format!(
                "Frozen {:?} is not supported",
                self.problem
            )));
        }
        let solver = self.solver()?;
        let mixtures = self.weights(&solver)?;
        let mut cases = vec![];
        let mut push = |weights: &[f64], input| {
            cases.push(Case {
                deck: self,
                index: cases.len() + 1,
                weights: weights.to_vec(),
                input,
            })
        };
        match self.problem {
            Problem::Rocket => {
                if self.temperatures_k.len() > 1 || self.enthalpies_j_kg.len() > 1 {
                    return Err(Error::Unsupported(
                        "Rocket t and h/r each require a single chamber value".into(),
                    ));
                }
                if !self.temperatures_k.is_empty() && !self.enthalpies_j_kg.is_empty() {
                    return Err(invalid("Rocket t and h/r are mutually exclusive"));
                }
                for frozen in [false, true] {
                    if (frozen && !self.frozen) || (!frozen && !self.equilibrium && self.frozen) {
                        continue;
                    }
                    for &pressure_pa in &self.pressures_pa {
                        for weights in &mixtures {
                            push(
                                weights,
                                CaseInput::Rocket {
                                    pressure_pa,
                                    freeze_station: if frozen {
                                        self.rocket.freeze_station
                                    } else {
                                        0
                                    },
                                },
                            );
                        }
                    }
                }
            }
            Problem::Shock => {
                if mixtures.len() != 1 {
                    return Err(Error::Unsupported(
                        "Shock mixture ratio schedules are not supported by NASA's deck interface"
                            .into(),
                    ));
                }
                let temps = if self.temperatures_k.is_empty() {
                    vec![
                        self.reactants[0]
                            .temperature_k
                            .ok_or_else(|| invalid("Shock requires t or a reactant temperature"))?,
                    ]
                } else {
                    self.temperatures_k.clone()
                };
                let models = if self.equilibrium && self.frozen {
                    vec![false, true]
                } else {
                    vec![!self.equilibrium]
                };
                for &incident_frozen in &models {
                    // NASA prints frozen reflection before equilibrium reflection.
                    let reflections = if self.reflected && models.len() == 2 {
                        vec![true, false]
                    } else {
                        vec![incident_frozen]
                    };
                    for reflected_frozen in reflections {
                        for (i, &speed) in self.shock_speeds.iter().enumerate() {
                            push(
                                &mixtures[0],
                                CaseInput::Shock {
                                    temperature_k: temps[i.min(temps.len() - 1)],
                                    pressure_pa: self.pressures_pa
                                        [i.min(self.pressures_pa.len() - 1)],
                                    speed,
                                    incident_frozen,
                                    reflected_frozen,
                                },
                            );
                        }
                    }
                }
            }
            Problem::Detonation => {
                for weights in &mixtures {
                    for &temperature_k in &self.temperatures_k {
                        for &pressure_pa in &self.pressures_pa {
                            push(
                                weights,
                                CaseInput::Detonation {
                                    temperature_k,
                                    pressure_pa,
                                },
                            );
                        }
                    }
                }
            }
            _ => {
                let state2 = if matches!(self.problem, Problem::TV | Problem::UV | Problem::SV) {
                    &self.volumes_m3_kg
                } else {
                    &self.pressures_pa
                };
                for weights in &mixtures {
                    let h;
                    let state1 = match self.problem {
                        Problem::TP | Problem::TV => &self.temperatures_k,
                        Problem::HP if self.enthalpies_j_kg.is_empty() => {
                            h = vec![solver.reactant_enthalpy(weights, &self.temperatures())?];
                            &h
                        }
                        Problem::HP => &self.enthalpies_j_kg,
                        Problem::UV => &self.energies_j_kg,
                        _ => &self.entropies_j_kg_k,
                    };
                    for &b in state2 {
                        for &a in state1 {
                            let mode = match self.problem {
                                Problem::TP => Mode::TP(a, b),
                                Problem::HP => Mode::HP(a, b),
                                Problem::SP => Mode::SP(a, b),
                                Problem::TV => Mode::TV(a, b),
                                Problem::UV => Mode::UV(a, b),
                                Problem::SV => Mode::SV(a, b),
                                _ => unreachable!(),
                            };
                            push(weights, CaseInput::Equilibrium(mode));
                        }
                    }
                }
            }
        }
        if cases.is_empty() {
            return Err(invalid("No cases in deck"));
        }
        Ok(cases)
    }
}

impl Case<'_> {
    pub fn run(&self) -> Result<DeckResult> {
        self.solve().map_err(|error| {
            let prefix = format!("case {:?} #{}: ", self.deck.name, self.index);
            match error {
                Error::Input(m) => Error::Input(prefix + &m),
                Error::Domain(m) => Error::Domain(prefix + &m),
                Error::Data(m) => Error::Data(prefix + &m),
                Error::Convergence(m) => Error::Convergence(prefix + &m),
                Error::Unsupported(m) => Error::Unsupported(prefix + &m),
            }
        })
    }

    fn solve(&self) -> Result<DeckResult> {
        let d = self.deck;
        let solver = d.solver()?;
        let output = match self.input {
            CaseInput::Equilibrium(mode) => {
                CaseOutput::Equilibrium(Box::new(solver.solve(mode, &self.weights)?))
            }
            CaseInput::Rocket {
                pressure_pa,
                freeze_station,
            } => {
                let r = &d.rocket;
                CaseOutput::Rocket(Box::new(solver.rocket_stations(
                    &self.weights,
                    &d.temperatures(),
                    pressure_pa,
                    &RocketOptions {
                        pi_p: &r.pressure_ratios,
                        subar: &r.subsonic_area_ratios,
                        supar: &r.supersonic_area_ratios,
                        iac: !r.finite_area,
                        ac_at: r.ac_at,
                        mdot: r.mass_flux_kg_s_m2,
                        n_frz: freeze_station,
                        hc: d.enthalpies_j_kg.first().copied(),
                        tc: d.temperatures_k.first().copied(),
                        tc_est: r.temperature_estimate_k,
                        transport: d.output.transport,
                    },
                )?))
            }
            CaseInput::Shock {
                temperature_k,
                pressure_pa,
                speed,
                incident_frozen,
                reflected_frozen,
            } => CaseOutput::Shock(Box::new(solver.shock(
                &self.weights,
                temperature_k,
                pressure_pa,
                speed,
                d.reflected,
                incident_frozen,
                reflected_frozen,
            )?)),
            CaseInput::Detonation {
                temperature_k,
                pressure_pa,
            } => CaseOutput::Detonation(Box::new(solver.detonation(
                &self.weights,
                temperature_k,
                pressure_pa,
                false,
            )?)),
        };
        let transport = if d.output.transport {
            match &output {
                CaseOutput::Equilibrium(st) => vec![st.transport()?],
                CaseOutput::Detonation(r) => vec![r.burned.transport()?],
                CaseOutput::Rocket(r) => r
                    .stations
                    .iter()
                    .map(|s| {
                        s.transport
                            .ok_or_else(|| Error::Data("Missing requested rocket transport".into()))
                    })
                    .collect::<Result<_>>()?,
                CaseOutput::Shock(r) => [&r.station1, &r.station2]
                    .into_iter()
                    .chain(r.station5.as_ref())
                    .map(|s| s.state.transport())
                    .collect::<Result<_>>()?,
            }
        } else {
            vec![]
        };
        let volume = d
            .reactants
            .iter()
            .zip(&self.weights)
            .map(|(r, w)| r.density_kg_m3.map(|rho| w / rho))
            .sum::<Option<f64>>();
        Ok(DeckResult {
            case_name: d.name.clone(),
            case_index: self.index,
            weights: self.weights.clone(),
            input: self.input,
            output,
            transport,
            reactant_density_kg_m3: volume.map(|v| self.weights.iter().sum::<f64>() / v),
        })
    }
}

/// Run every case in a text containing one or more decks. Errors identify the
/// failing case; no partial list is returned. Use `parse_deck`/`Deck::cases`
/// to collect successes and failures independently.
pub fn run_deck(text: &str) -> Result<Vec<DeckResult>> {
    let mut results = vec![];
    for deck in parse_deck(text)? {
        for case in deck.cases()? {
            results.push(case.run()?);
        }
    }
    Ok(results)
}
