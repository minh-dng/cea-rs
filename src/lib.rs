//! Experimental native Rust chemical equilibrium. No Fortran runtime is used.
#![allow(clippy::needless_range_loop)]
#[allow(dead_code)]
mod equilibrium_core;
#[cfg(feature = "python")]
mod python;
pub mod rocket;
pub mod solver;
pub mod thermo;

/// NASA CEA gas constant, J/(kmol K), kept for oracle comparisons.
pub const R: f64 = 8314.51;
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Clone)]
pub enum Error {
    Input(String),
    Domain(String),
    Data(String),
    Convergence(String),
    Unsupported(String),
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (kind, message) = match self {
            Self::Input(m) => ("invalid_input", m),
            Self::Domain(m) => ("out_of_domain", m),
            Self::Data(m) => ("invalid_data", m),
            Self::Convergence(m) => ("not_converged", m),
            Self::Unsupported(m) => ("unsupported", m),
        };
        write!(f, "{kind}: {message}")
    }
}
impl std::error::Error for Error {}
