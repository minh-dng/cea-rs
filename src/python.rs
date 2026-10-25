//! Python boundary: SI units, typed failures and a distinct import namespace.
use crate::{
    Error,
    rocket::{FlowMode, RocketResult, Station},
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
#[pymethods]
impl Reactant {
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
    #[new]
    #[pyo3(signature=(reactants, *, condensed=true, ions=false, only=None, insert=Vec::new()))]
    fn new(
        py: Python<'_>,
        reactants: Vec<Py<Reactant>>,
        condensed: bool,
        ions: bool,
        only: Option<Vec<String>>,
        insert: Vec<String>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: Solver::new(
                reactants.iter().map(|r| r.borrow(py).clone()).collect(),
                condensed,
                ions,
                only,
                insert,
            )?,
        })
    }
    /// Solve TP/HP/SP/TV/UV/SV in SI (K, Pa, J/kg, J/kg/K, m³/kg).
    fn solve(
        &self,
        py: Python<'_>,
        mode: &str,
        state1: f64,
        state2: f64,
        weights: Vec<f64>,
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
        Ok(py.detach(|| self.inner.solve(mode, &weights))?)
    }
    /// Return stream enthalpy in J/kg at specified reactant temperatures.
    fn reactant_enthalpy(&self, weights: Vec<f64>, temperatures_k: Vec<f64>) -> PyResult<f64> {
        Ok(self.inner.reactant_enthalpy(&weights, &temperatures_k)?)
    }
    #[pyo3(signature=(weights,temperatures_k,chamber_pressure_pa,expansion_ratio, *, flow="equilibrium",iac=true,transport=false))]
    fn rocket(
        &self,
        py: Python<'_>,
        weights: Vec<f64>,
        temperatures_k: Vec<f64>,
        chamber_pressure_pa: f64,
        expansion_ratio: f64,
        flow: &str,
        iac: bool,
        transport: bool,
    ) -> PyResult<RocketResult> {
        if !iac || transport {
            return Err(PyNotImplementedError::new_err(
                "Finite-area chambers and transport are not implemented",
            ));
        }
        let flow = match flow {
            "equilibrium" => FlowMode::Equilibrium,
            "frozen_at_chamber" => FlowMode::FrozenAtChamber,
            "frozen_at_throat" => FlowMode::FrozenAtThroat,
            _ => return Err(PyValueError::new_err("Unknown flow mode")),
        };
        Ok(py.detach(|| {
            self.inner.rocket(
                &weights,
                &temperatures_k,
                chamber_pressure_pa,
                expansion_ratio,
                flow,
            )
        })?)
    }
}
#[pymethods]
impl State {
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
    fn cp_eq(&self) -> Option<f64> {
        self.equilibrium.then_some(self.cp_equilibrium_j_kg_k)
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
}
#[pymethods]
impl RocketResult {
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
    m.add_class::<Reactant>()?;
    m.add_class::<PySolver>()?;
    m.add_class::<State>()?;
    m.add_class::<RocketResult>()?;
    m.add_class::<Station>()?;
    m.add("ConvergenceError", m.py().get_type::<ConvergenceError>())?;
    m.add("DomainError", m.py().get_type::<DomainError>())?;
    m.add("__version__", "0.1.0a1")?;
    m.add("R", crate::R)?;
    Ok(())
}
