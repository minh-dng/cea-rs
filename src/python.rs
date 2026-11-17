//! Python boundary: SI units, typed failures and a distinct import namespace.
use crate::{
    Error,
    deck::{CaseOutput, DeckResult},
    detonation::{DetonationResult, UnburnedState},
    rocket::{RocketOptions, RocketResult, Station},
    solver::{Mode, Solver, State},
    thermo::Reactant,
};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::collections::BTreeMap;
pyo3::create_exception!(_cea_rs, ConvergenceError, PyRuntimeError);
pyo3::create_exception!(_cea_rs, DomainError, PyValueError);
impl From<Error> for PyErr {
    fn from(e: Error) -> Self {
        match e {
            Error::Input(_) | Error::Data(_) => PyValueError::new_err(e.to_string()),
            Error::Domain(_) => DomainError::new_err(e.to_string()),
            Error::Convergence(_) => ConvergenceError::new_err(e.to_string()),
            Error::Unsupported(_) => PyNotImplementedError::new_err(e.to_string()),
        }
    }
}

#[pyfunction]
fn run_deck(py: Python<'_>, text: &str) -> PyResult<Vec<DeckResult>> {
    Ok(py.detach(|| crate::deck::run_deck(text))?)
}

type NasaSpecies = (String, f64, i32, BTreeMap<String, f64>);

/// Read-only database facts in NASA's compiled gas-then-condensed name order.
#[pyfunction]
fn _nasa_species() -> Vec<NasaSpecies> {
    let mut species = thermo_species::Species::all().to_vec();
    species.sort_by_key(|sp| (sp.data().phase() != 0, sp.data().symbol()));
    species
        .iter()
        .map(|sp| {
            let d = sp.data();
            (
                d.symbol().into(),
                d.mw(),
                d.phase() as i32,
                d.constituents()
                    .iter()
                    .map(|(n, e)| (format!("{e:?}"), *n))
                    .collect(),
            )
        })
        .collect()
}

#[pymethods]
impl DeckResult {
    #[getter]
    fn case_name(&self) -> &str {
        &self.case_name
    }
    #[getter]
    fn case_index(&self) -> usize {
        self.case_index
    }
    #[getter]
    fn weights(&self) -> Vec<f64> {
        self.weights.clone()
    }
    #[getter]
    fn kind(&self) -> &'static str {
        match self.output {
            CaseOutput::Equilibrium(_) => "equilibrium",
            CaseOutput::Rocket(_) => "rocket",
            CaseOutput::Shock(_) => "shock",
            CaseOutput::Detonation(_) => "detonation",
        }
    }
    #[getter]
    fn result(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(match &self.output {
            CaseOutput::Equilibrium(r) => Py::new(py, r.as_ref().clone())?.into_any(),
            CaseOutput::Rocket(r) => Py::new(py, r.as_ref().clone())?.into_any(),
            CaseOutput::Shock(r) => Py::new(py, r.as_ref().clone())?.into_any(),
            CaseOutput::Detonation(r) => Py::new(py, r.as_ref().clone())?.into_any(),
        })
    }
    #[getter]
    fn transport(&self) -> Vec<crate::transport::Transport> {
        self.transport.clone()
    }
    #[getter]
    fn reactant_density_kg_m3(&self) -> Option<f64> {
        self.reactant_density_kg_m3
    }
}
#[pymethods]
impl Reactant {
    #[getter]
    fn _nasa_formula(&self) -> BTreeMap<String, f64> {
        self.formula
            .iter()
            .map(|(n, e)| (format!("{e:?}"), *n))
            .collect()
    }
    #[getter]
    fn _nasa_phase(&self) -> i32 {
        self.phase.into()
    }
    #[pyo3(name = "get_valid_temperature_range")]
    fn py_valid_temperature_range(&self) -> (f64, f64) {
        self.get_valid_temperature_range()
    }
    #[new]
    #[pyo3(signature=(name, *, formula=None, molecular_weight=None, enthalpy_j_mol=None, reference_temperature_k=298.15))]
    fn py_new(
        name: String,
        formula: Option<BTreeMap<String, f64>>,
        molecular_weight: Option<f64>,
        enthalpy_j_mol: Option<f64>,
        reference_temperature_k: f64,
    ) -> PyResult<Self> {
        if let Some(formula) = formula {
            Ok(Self::custom(
                name,
                formula,
                molecular_weight.ok_or_else(|| {
                    PyValueError::new_err("Custom reactants require molecular_weight")
                })?,
                enthalpy_j_mol.ok_or_else(|| {
                    PyValueError::new_err("Custom reactants require enthalpy_j_mol")
                })?,
                reference_temperature_k,
            )?)
        } else {
            if molecular_weight.is_some()
                || enthalpy_j_mol.is_some()
                || reference_temperature_k != 298.15
            {
                return Err(PyValueError::new_err(
                    "Named NASA reactants cannot override reference data; supply a formula for custom data",
                ));
            }
            Ok(Self::named(&name)?)
        }
    }
    #[getter(name)]
    fn get_name(&self) -> String {
        self.name.clone()
    }
    #[getter(molecular_weight)]
    fn get_molecular_weight(&self) -> f64 {
        self.molecular_weight
    }
    /// Return Cp (J/mol/K), H (J/mol), S° (J/mol/K); unavailable values are None.
    fn thermo(&self, temperature_k: f64) -> PyResult<(Option<f64>, f64, Option<f64>)> {
        Ok(self.properties(temperature_k)?)
    }
}

#[pyclass(name = "Solver")]
struct PySolver {
    inner: Solver,
}
#[pymethods]
impl PySolver {
    fn of_ratio_to_weights(
        &self,
        oxidant_weights: Vec<f64>,
        fuel_weights: Vec<f64>,
        of_ratio: f64,
    ) -> PyResult<Vec<f64>> {
        Ok(self
            .inner
            .of_ratio_to_weights(&oxidant_weights, &fuel_weights, of_ratio)?)
    }
    fn chem_eq_ratio_to_of_ratio(
        &self,
        oxidant_weights: Vec<f64>,
        fuel_weights: Vec<f64>,
        chem_eq_ratio: f64,
    ) -> PyResult<f64> {
        Ok(self
            .inner
            .chem_eq_ratio_to_of_ratio(&oxidant_weights, &fuel_weights, chem_eq_ratio)?)
    }
    fn of_ratio_to_chem_eq_ratio(
        &self,
        oxidant_weights: Vec<f64>,
        fuel_weights: Vec<f64>,
        of_ratio: f64,
    ) -> PyResult<f64> {
        Ok(self
            .inner
            .of_ratio_to_chem_eq_ratio(&oxidant_weights, &fuel_weights, of_ratio)?)
    }
    fn weight_eq_ratio_to_of_ratio(
        &self,
        oxidant_weights: Vec<f64>,
        fuel_weights: Vec<f64>,
        weight_eq_ratio: f64,
    ) -> PyResult<f64> {
        Ok(self.inner.weight_eq_ratio_to_of_ratio(
            &oxidant_weights,
            &fuel_weights,
            weight_eq_ratio,
        )?)
    }
    fn of_ratio_to_weight_eq_ratio(
        &self,
        oxidant_weights: Vec<f64>,
        fuel_weights: Vec<f64>,
        of_ratio: f64,
    ) -> PyResult<f64> {
        Ok(self
            .inner
            .of_ratio_to_weight_eq_ratio(&oxidant_weights, &fuel_weights, of_ratio)?)
    }
    fn weights_to_moles(&self, weights: Vec<f64>) -> PyResult<Vec<f64>> {
        Ok(self.inner.weights_to_moles(&weights)?)
    }
    fn moles_to_weights(&self, moles: Vec<f64>) -> PyResult<Vec<f64>> {
        Ok(self.inner.moles_to_weights(&moles)?)
    }
    fn get_valence(&self) -> PyResult<BTreeMap<String, (f64, f64)>> {
        Ok(self.inner.get_valence()?)
    }
    #[pyo3(signature=(mode,state1,state2,weights, *, finite_difference=false,relative_step=1e-4))]
    #[allow(clippy::too_many_arguments)]
    fn derivatives(
        &self,
        py: Python<'_>,
        mode: &str,
        state1: f64,
        state2: f64,
        weights: Vec<f64>,
        finite_difference: bool,
        relative_step: f64,
    ) -> PyResult<crate::derivatives::EqDerivatives> {
        let mode = match mode {
            "TP" => Mode::TP(state1, state2),
            "HP" => Mode::HP(state1, state2),
            "SP" => Mode::SP(state1, state2),
            "TV" => Mode::TV(state1, state2),
            "UV" => Mode::UV(state1, state2),
            "SV" => Mode::SV(state1, state2),
            _ => {
                return Err(PyValueError::new_err(
                    "Supported equilibrium modes: TP HP SP TV UV SV",
                ));
            }
        };
        Ok(py.detach(|| {
            if finite_difference {
                self.inner
                    .finite_difference_derivatives(mode, &weights, relative_step)
            } else {
                self.inner.derivatives(mode, &weights)
            }
        })?)
    }
    #[new]
    #[pyo3(signature=(reactants, *, condensed=true, ions=false, only=None, insert=Vec::new(), omit=Vec::new(), trace=0.0))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        reactants: Vec<Py<Reactant>>,
        condensed: bool,
        ions: bool,
        only: Option<Vec<String>>,
        insert: Vec<String>,
        omit: Vec<String>,
        trace: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: Solver::new(
                reactants.iter().map(|r| r.borrow(py).clone()).collect(),
                condensed,
                ions,
                only,
                insert,
            )?
            .with_omit(omit)?
            .with_trace(trace)?,
        })
    }
    /// Solve TP/HP/SP/TV/UV/SV in SI (K, Pa, J/kg, J/kg/K, m³/kg).
    #[pyo3(signature=(mode, state1, state2, weights, *, start=None, transport=false))]
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &self,
        py: Python<'_>,
        mode: &str,
        state1: f64,
        state2: f64,
        weights: Vec<f64>,
        start: Option<PyRef<'_, State>>,
        transport: bool,
    ) -> PyResult<State> {
        let mode = match mode {
            "TP" => Mode::TP(state1, state2),
            "HP" => Mode::HP(state1, state2),
            "SP" => Mode::SP(state1, state2),
            "TV" => Mode::TV(state1, state2),
            "UV" => Mode::UV(state1, state2),
            "SV" => Mode::SV(state1, state2),
            _ => {
                return Err(PyNotImplementedError::new_err(
                    "Supported equilibrium modes: TP HP SP TV UV SV",
                ));
            }
        };
        let start = start.map(|s| s.clone());
        let mut solver = self.inner.clone();
        solver.transport = transport;
        Ok(py.detach(|| {
            solver.solve_with_start(
                mode,
                &weights,
                start.as_ref().and_then(|s| s.equilibrium_seed.as_ref()),
            )
        })?)
    }
    /// Return stream enthalpy in J/kg at specified reactant temperatures.
    fn reactant_enthalpy(&self, weights: Vec<f64>, temperatures_k: Vec<f64>) -> PyResult<f64> {
        Ok(self.inner.reactant_enthalpy(&weights, &temperatures_k)?)
    }
    /// Solve a CJ detonation from mass weights, unburned K and Pa.
    #[pyo3(signature=(weights, temperature_k, pressure_pa, *, frozen=false))]
    fn detonation(
        &self,
        py: Python<'_>,
        weights: Vec<f64>,
        temperature_k: f64,
        pressure_pa: f64,
        frozen: bool,
    ) -> PyResult<DetonationResult> {
        Ok(py.detach(|| {
            self.inner
                .detonation(&weights, temperature_k, pressure_pa, frozen)
        })?)
    }
    /// Incident/reflected shock: mass weights, K, Pa, and exactly one speed input.
    #[pyo3(signature=(weights,temperature_k,pressure_pa, *, u1=None,mach1=None,reflected=false,incident_frozen=false,reflected_frozen=false,transport=false))]
    #[allow(clippy::too_many_arguments)]
    fn shock(
        &self,
        py: Python<'_>,
        weights: Vec<f64>,
        temperature_k: f64,
        pressure_pa: f64,
        u1: Option<f64>,
        mach1: Option<f64>,
        reflected: bool,
        incident_frozen: bool,
        reflected_frozen: bool,
        transport: bool,
    ) -> PyResult<crate::shock::ShockResult> {
        use crate::shock::ShockSpeed;
        let speed = match (u1, mach1) {
            (Some(u), None) => ShockSpeed::Velocity(u),
            (None, Some(m)) => ShockSpeed::Mach(m),
            _ => return Err(PyValueError::new_err("Specify exactly one of u1 or mach1")),
        };
        let mut solver = self.inner.clone();
        solver.transport = transport;
        Ok(py.detach(|| {
            solver.shock(
                &weights,
                temperature_k,
                pressure_pa,
                speed,
                reflected,
                incident_frozen,
                reflected_frozen,
            )
        })?)
    }
    #[pyo3(signature=(weights,temperatures_k,chamber_pressure_pa,expansion_ratio=None, *, flow="equilibrium",iac=true,transport=false,pi_p=Vec::new(),subar=Vec::new(),supar=Vec::new(),ac_at=None,mdot=None,n_frz=None,hc=None,tc=None,tc_est=None))]
    #[allow(clippy::too_many_arguments)]
    fn rocket(
        &self,
        py: Python<'_>,
        weights: Vec<f64>,
        temperatures_k: Vec<f64>,
        chamber_pressure_pa: f64,
        expansion_ratio: Option<f64>,
        flow: &str,
        iac: bool,
        transport: bool,
        pi_p: Vec<f64>,
        subar: Vec<f64>,
        mut supar: Vec<f64>,
        ac_at: Option<f64>,
        mdot: Option<f64>,
        n_frz: Option<isize>,
        hc: Option<f64>,
        tc: Option<f64>,
        tc_est: Option<f64>,
    ) -> PyResult<RocketResult> {
        if !iac && ac_at.is_none() && mdot.is_none() && expansion_ratio.is_some() {
            return Err(PyNotImplementedError::new_err(
                "FAC calls require ac_at or mdot",
            ));
        }
        let n_frz = n_frz
            .map(|n| {
                usize::try_from(n).map_err(|_| {
                    PyValueError::new_err("n_frz must be zero or a positive station index")
                })
            })
            .transpose()?;
        let freeze = match flow {
            "equilibrium" => 0,
            "frozen_at_chamber" => {
                if iac {
                    1
                } else {
                    3
                }
            }
            "frozen_at_throat" => {
                if iac {
                    2
                } else {
                    4
                }
            }
            _ => return Err(PyValueError::new_err("Unknown flow mode")),
        };
        if n_frz.is_some_and(|n| freeze != 0 && n != freeze) {
            return Err(PyValueError::new_err("flow conflicts with n_frz"));
        }
        if let Some(area) = expansion_ratio {
            if !supar.is_empty() {
                return Err(PyValueError::new_err(
                    "Use expansion_ratio or supar, not both",
                ));
            }
            supar.push(area);
        }
        let result = py.detach(|| {
            self.inner.rocket_stations(
                &weights,
                &temperatures_k,
                chamber_pressure_pa,
                &RocketOptions {
                    pi_p: &pi_p,
                    subar: &subar,
                    supar: &supar,
                    iac,
                    ac_at,
                    mdot,
                    n_frz: n_frz.unwrap_or(freeze),
                    hc,
                    tc,
                    tc_est,
                    transport,
                },
            )
        })?;
        if expansion_ratio.is_some() && !result.complete {
            return Err(DomainError::new_err(
                result
                    .stop_reason
                    .unwrap_or_else(|| "Frozen exit was not completed".into()),
            ));
        }
        Ok(result)
    }
}
#[pymethods]
impl State {
    #[getter]
    fn converged(&self) -> bool {
        self.converged
    }
    #[getter]
    fn last_error(&self) -> i32 {
        self.last_error
    }
    #[getter]
    fn _nasa_n(&self) -> f64 {
        self.gas_moles
    }
    #[getter]
    fn _nasa_post_processed(&self) -> bool {
        self.post_processed
    }
    #[getter]
    fn _nasa_nj(&self) -> BTreeMap<String, f64> {
        self.amounts
            .iter()
            .map(|(n, sp)| (sp.data().symbol().into(), *n))
            .collect()
    }
    #[getter]
    fn _nasa_ln_nj(&self) -> Option<BTreeMap<String, f64>> {
        self.equilibrium_seed.as_ref().map(|seed| {
            assert_eq!(seed.species.len(), seed.ln_nj.len());
            seed.species
                .iter()
                .zip(&seed.ln_nj)
                .filter(|(sp, _)| sp.data().phase() == 0)
                .map(|(sp, ln_n)| (sp.data().symbol().into(), *ln_n))
                .collect()
        })
    }
    #[getter(gibbs_energy_j_kg)]
    fn gibbs_energy_py(&self) -> f64 {
        crate::mixture::gibbs_energy(self)
    }
    #[getter(cv_frozen_j_kg_k)]
    fn cv_frozen_py(&self) -> f64 {
        crate::mixture::cv_frozen(self)
    }
    #[getter(cv_equilibrium_j_kg_k)]
    fn cv_equilibrium_py(&self) -> f64 {
        self.cv_equilibrium_j_kg_k
    }
    #[getter(specific_volume_m3_kg)]
    fn specific_volume_py(&self) -> f64 {
        crate::mixture::specific_volume(self)
    }
    #[getter(M)]
    fn nasa_m_py(&self) -> f64 {
        crate::mixture::molecular_weight_m(self)
    }
    #[getter(MW)]
    fn nasa_mw_py(&self) -> f64 {
        crate::mixture::molecular_weight_mw(self)
    }
    #[pyo3(name = "transport")]
    fn transport_py(&self, py: Python<'_>) -> PyResult<crate::transport::Transport> {
        Ok(py.detach(|| self.transport())?)
    }
    #[getter(temperature_k)]
    fn t(&self) -> f64 {
        self.temperature_k
    }
    #[getter(pressure_pa)]
    fn p(&self) -> f64 {
        self.pressure_pa
    }
    #[getter(density_kg_m3)]
    fn rho(&self) -> f64 {
        self.density_kg_m3
    }
    #[getter(enthalpy_j_kg)]
    fn h(&self) -> f64 {
        self.enthalpy_j_kg
    }
    #[getter(energy_j_kg)]
    fn u(&self) -> f64 {
        self.energy_j_kg
    }
    #[getter(entropy_j_kg_k)]
    fn s(&self) -> f64 {
        self.entropy_j_kg_k
    }
    #[getter(cp_equilibrium_j_kg_k)]
    fn cp_eq(&self) -> f64 {
        self.cp_equilibrium_j_kg_k
    }
    #[getter(cp_frozen_j_kg_k)]
    fn cp_fr(&self) -> f64 {
        self.cp_frozen_j_kg_k
    }
    #[getter(gamma_s)]
    fn gamma(&self) -> f64 {
        self.gamma_s
    }
    #[getter(molecular_weight_kg_kmol)]
    fn mw(&self) -> f64 {
        self.molecular_weight_kg_kmol
    }
    #[getter(equilibrium)]
    fn equilibrium_py(&self) -> bool {
        self.equilibrium
    }
    #[getter(extrapolated_species)]
    fn get_extrapolated_species(&self) -> Vec<String> {
        self.extrapolated_species.clone()
    }
    #[getter(trace_extrapolated_species)]
    fn trace_extrapolated(&self) -> Vec<String> {
        self.trace_extrapolated_species.clone()
    }
    #[getter(iterations)]
    fn iters(&self) -> usize {
        self.iterations
    }
    #[getter(element_residual)]
    fn elements(&self) -> f64 {
        self.element_residual
    }
    #[getter(constraint_residual)]
    fn constraint(&self) -> f64 {
        self.constraint_residual
    }
    #[getter(mole_fractions)]
    fn x(&self) -> BTreeMap<String, f64> {
        self.mole_fractions.clone()
    }
    #[getter(mass_fractions)]
    fn y(&self) -> BTreeMap<String, f64> {
        self.mass_fractions.clone()
    }
    #[getter(sound_speed_m_s)]
    fn a(&self) -> f64 {
        self.sound_speed()
    }
}
#[pymethods]
impl Station {
    #[getter(name)]
    fn name_py(&self) -> &str {
        self.name
    }
    #[getter(sound_speed_m_s)]
    fn sonic(&self) -> f64 {
        self.sound_speed_m_s
    }
    #[getter(cstar_m_s)]
    fn cstar_py(&self) -> f64 {
        self.cstar_m_s
    }
    #[getter(thrust_coefficient)]
    fn cf(&self) -> f64 {
        self.thrust_coefficient
    }
    #[getter(isp_m_s)]
    fn isp(&self) -> f64 {
        self.isp_m_s
    }
    #[getter(isp_vacuum_m_s)]
    fn ivac(&self) -> f64 {
        self.isp_vacuum_m_s
    }
    #[getter(state)]
    fn state_py(&self) -> State {
        self.state.clone()
    }
    #[getter(velocity_m_s)]
    fn velocity(&self) -> f64 {
        self.velocity_m_s
    }
    #[getter(mach)]
    fn mach_py(&self) -> f64 {
        self.mach
    }
    #[getter(area_ratio)]
    fn area(&self) -> Option<f64> {
        self.area_ratio
    }
    #[getter(transport)]
    fn transport_py(&self) -> Option<crate::transport::Transport> {
        self.transport
    }
}
#[pymethods]
impl RocketResult {
    #[getter(stations)]
    fn stations_py(&self) -> Vec<Station> {
        self.stations.clone()
    }
    #[getter(throat_index)]
    fn throat_index_py(&self) -> usize {
        self.throat_index
    }
    #[getter(complete)]
    fn complete_py(&self) -> bool {
        self.complete
    }
    #[getter(stop_reason)]
    fn stop_reason_py(&self) -> Option<String> {
        self.stop_reason.clone()
    }
    #[getter(chamber)]
    fn chamber_py(&self) -> Station {
        self.chamber.clone()
    }
    #[getter(throat)]
    fn throat_py(&self) -> Station {
        self.throat.clone()
    }
    #[getter(exit)]
    fn exit_py(&self) -> Station {
        self.exit.clone()
    }
    #[getter(cstar_m_s)]
    fn cstar(&self) -> f64 {
        self.cstar_m_s
    }
    #[getter(vacuum_velocity_m_s)]
    fn vacuum(&self) -> f64 {
        self.vacuum_velocity_m_s
    }
}
#[pymodule]
fn _cea_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_nasa_species, m)?)?;
    m.add_class::<DeckResult>()?;
    m.add_function(wrap_pyfunction!(run_deck, m)?)?;
    m.add_class::<crate::derivatives::EqDerivatives>()?;
    m.add_class::<Reactant>()?;
    m.add_class::<PySolver>()?;
    m.add_class::<State>()?;
    m.add_class::<RocketResult>()?;
    m.add_class::<Station>()?;
    m.add_class::<crate::transport::Transport>()?;
    m.add_class::<DetonationResult>()?;
    m.add_class::<UnburnedState>()?;
    m.add_class::<crate::shock::ShockStation>()?;
    m.add_class::<crate::shock::ShockResult>()?;
    m.add("ConvergenceError", m.py().get_type::<ConvergenceError>())?;
    m.add("DomainError", m.py().get_type::<DomainError>())?;
    m.add("__version__", "0.1.0a3")?;
    m.add("R", crate::R)?;
    Ok(())
}

#[pymethods]
impl crate::derivatives::EqDerivatives {
    #[getter(temperature)]
    fn temperature_py(&self) -> Vec<f64> {
        self.temperature.clone()
    }
    #[getter(total_moles)]
    fn total_moles_py(&self) -> Vec<f64> {
        self.total_moles.clone()
    }
    #[getter(species_amounts)]
    fn species_amounts_py(&self) -> BTreeMap<String, Vec<f64>> {
        self.species_amounts.clone()
    }
    #[getter(enthalpy)]
    fn enthalpy_py(&self) -> Vec<f64> {
        self.enthalpy.clone()
    }
    #[getter(energy)]
    fn energy_py(&self) -> Vec<f64> {
        self.energy.clone()
    }
    #[getter(gibbs_energy)]
    fn gibbs_energy_py(&self) -> Vec<f64> {
        self.gibbs_energy.clone()
    }
    #[getter(entropy)]
    fn entropy_py(&self) -> Vec<f64> {
        self.entropy.clone()
    }
    #[getter(closure_defect)]
    fn closure_defect_py(&self) -> Option<f64> {
        self.closure_defect
    }
}
