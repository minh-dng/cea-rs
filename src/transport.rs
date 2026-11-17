//! NASA CEA gas transport, derived from `transport.f90` and `equilibrium.f90`
//! at 2f79a647e85737742fb08a127da00d1ea1d069da, under Apache-2.0.
//!
//! Condensed species do not participate. NASA's end fits are extrapolated;
//! missing data use its kinetic-theory and rigid-sphere estimates. Reaction
//! transport uses a gas component basis, not the bulk equilibrium heat capacity.
use crate::{Error, R, Result, solver::State};
use nalgebra::{DMatrix, DVector};
use std::{collections::BTreeSet, f64::consts::PI, sync::OnceLock};
use thermo_species::{Constituent, Species};

/// Gas transport properties in SI, including NASA's reaction contribution.
/// The heat capacities are NASA's transport values: gas species only, over the
/// transport species selection, so they can differ from the state's.
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object)
)]
#[derive(Debug, Clone, Copy, Default)]
pub struct Transport {
    pub viscosity_pa_s: f64,
    pub cp_frozen_j_kg_k: f64,
    pub cp_equilibrium_j_kg_k: f64,
    pub conductivity_frozen_w_m_k: f64,
    pub conductivity_equilibrium_w_m_k: f64,
    pub prandtl_frozen: f64,
    pub prandtl_equilibrium: f64,
}

#[derive(Debug)]
struct Fit {
    hi: f64,
    a: [f64; 4],
}

#[derive(Debug)]
struct Record {
    names: [String; 2],
    viscosity: Vec<Fit>,
    conductivity: Vec<Fit>,
}

fn field(line: &str, lo: usize, hi: usize) -> Result<&str> {
    line.get(lo..hi)
        .map(str::trim)
        .ok_or_else(|| Error::Data("Truncated NASA transport record".into()))
}

fn number(line: &str, lo: usize, hi: usize) -> Result<f64> {
    let value = field(line, lo, hi)?
        .replace(' ', "")
        .replace(['D', 'd'], "E")
        .parse::<f64>()
        .map_err(|_| Error::Data("Invalid NASA transport number".into()))?;
    if !value.is_finite() {
        return Err(Error::Data("Nonfinite NASA transport number".into()));
    }
    Ok(value)
}

fn parse(text: &str) -> Result<Vec<Record>> {
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some("transport property coefficients") {
        return Err(Error::Data("Missing NASA transport header".into()));
    }
    let mut records = Vec::new();
    let mut names_seen = BTreeSet::new();
    while let Some(header) = lines.next() {
        if header.trim() == "end" {
            if lines.any(|l| !l.trim().is_empty()) {
                return Err(Error::Data("Trailing NASA transport records".into()));
            }
            return Ok(records);
        }
        let names: [String; 2] = [field(header, 0, 16)?.into(), field(header, 16, 32)?.into()];
        if names[0].is_empty()
            || !names_seen.insert(names.clone())
            || field(header, 34, 35)? != "V"
            || field(header, 36, 37)? != "C"
        {
            return Err(Error::Data("Invalid NASA transport species header".into()));
        }
        let count = |lo| -> Result<usize> {
            let n = field(header, lo, lo + 1)?
                .parse::<usize>()
                .map_err(|_| Error::Data("Invalid NASA transport interval count".into()))?;
            if n > 3 {
                return Err(Error::Data("Too many NASA transport intervals".into()));
            }
            Ok(n)
        };
        let (nv, nc) = (count(35)?, count(37)?);
        let mut record = Record {
            names,
            viscosity: vec![],
            conductivity: vec![],
        };
        for _ in 0..nv + nc {
            let line = lines
                .next()
                .ok_or_else(|| Error::Data("Missing NASA transport fit".into()))?;
            let fits = match field(line, 1, 2)? {
                "V" => &mut record.viscosity,
                "C" => &mut record.conductivity,
                _ => return Err(Error::Data("Invalid NASA transport fit kind".into())),
            };
            let lo = number(line, 2, 11)?;
            let hi = number(line, 11, 20)?;
            if lo <= 0.0 || hi <= lo || fits.last().is_some_and(|f: &Fit| f.hi != lo) {
                return Err(Error::Data(
                    "Invalid NASA transport temperature interval".into(),
                ));
            }
            fits.push(Fit {
                hi,
                a: [
                    number(line, 20, 35)?,
                    number(line, 35, 50)?,
                    number(line, 50, 65)?,
                    number(line, 65, 80)?,
                ],
            });
        }
        if record.viscosity.len() != nv || record.conductivity.len() != nc {
            return Err(Error::Data("NASA transport interval count mismatch".into()));
        }
        records.push(record);
    }
    Err(Error::Data("Missing NASA transport end record".into()))
}

fn database() -> Result<&'static [Record]> {
    static DB: OnceLock<Result<Vec<Record>>> = OnceLock::new();
    DB.get_or_init(|| parse(include_str!("../vendor/nasa-transport/trans.inp")))
        .as_deref()
        .map_err(Clone::clone)
}

fn evaluate(fits: &[Fit], t: f64) -> f64 {
    // NASA returns log(value)=0 for an absent interval, before exponentiation.
    let Some(fit) = fits.iter().find(|f| t <= f.hi).or_else(|| fits.last()) else {
        return 1.0;
    };
    let [a, b, c, d] = fit.a;
    (a * t.ln() + b / t + c / (t * t) + d).exp()
}

fn coefficient(sp: Species, element: Constituent) -> f64 {
    sp.data()
        .constituents()
        .iter()
        .filter(|(_, e)| *e == element)
        .map(|(n, _)| n)
        .sum()
}

/// Rebuild NASA's abundance-selected, row-reduced gas component basis.
fn basis(gas: &[(f64, Species)], atomic: bool) -> Result<(Vec<usize>, DMatrix<f64>)> {
    let elements: BTreeSet<_> = gas
        .iter()
        .flat_map(|(_, s)| s.data().constituents())
        .map(|(_, e)| *e)
        .filter(|e| *e != Constituent::E)
        .collect();
    let elements: Vec<_> = elements.into_iter().collect();
    let ne = elements.len();
    let ng = gas.len();
    let mut a = DMatrix::from_fn(ne, ng, |i, j| coefficient(gas[j].1, elements[i]));
    let mut components = Vec::new();
    for i in 0..ne {
        let atomic = (0..ng).find(|&j| {
            (a[(i, j)] - 1.0).abs() < 1e-8
                && (gas[j]
                    .1
                    .data()
                    .constituents()
                    .iter()
                    .map(|(n, _)| n.abs())
                    .sum::<f64>()
                    - 1.0)
                    .abs()
                    < 1e-8
        });
        components.push(
            atomic
                .or_else(|| (0..ng).find(|&j| a[(i, j)].abs() > 1e-10))
                .ok_or_else(|| Error::Convergence("Missing gas transport component".into()))?,
        );
    }
    if atomic {
        return Ok((components, a));
    }
    let initial = components.clone();
    let mut work: Vec<_> = gas.iter().map(|(n, _)| *n).collect();
    let mut rows = Vec::new();
    while rows.len() < ne {
        let Some(j) = (0..ng).max_by(|&i, &j| work[i].total_cmp(&work[j]).then(j.cmp(&i))) else {
            break;
        };
        if work[j] <= 0.0 {
            break;
        }
        work[j] = -1.0;
        if coefficient(gas[j].1, Constituent::E).abs() > 1e-10 {
            continue;
        }
        for i in 0..ne {
            if a[(i, j)] <= 1e-10 || rows.contains(&i) {
                continue;
            }
            if rows
                .iter()
                .any(|&l| (0..ne).all(|r| a[(r, j)] == a[(r, components[l])]))
            {
                continue;
            }
            if (0..ne).any(|r| {
                r != i
                    && (a[(i, j)] * a[(r, initial[r])] - a[(i, initial[r])] * a[(r, j)]).abs()
                        <= 1e-10
            }) {
                continue;
            }
            components[i] = j;
            rows.push(i);
            break;
        }
    }
    // Normalize even an unchanged fallback basis (explicit molecular-only pools).
    let mut independent = Vec::new();
    for i in 0..ne {
        let mut j = components[i];
        if a[(i, j)] == 0.0 {
            // Dependent element rows do not add a gas component or a reaction.
            let Some(k) = (0..ng).find(|&k| a[(i, k)].abs() > 1e-10) else {
                continue;
            };
            j = k;
            components[i] = j;
        }
        let pivot = a[(i, j)];
        independent.push(i);
        for k in 0..ng {
            a[(i, k)] /= pivot;
        }
        for r in 0..ne {
            if r == i {
                continue;
            }
            let factor = a[(r, j)];
            for k in 0..ng {
                a[(r, k)] -= factor * a[(i, k)];
                if a[(r, k)].abs() < 1e-5 {
                    a[(r, k)] = 0.0;
                }
            }
        }
    }
    Ok((
        independent.iter().map(|&i| components[i]).collect(),
        a.select_rows(&independent),
    ))
}

/// Compute NASA's gas transport model from the state's composition.
/// As in NASA's frozen rocket calculations, the equilibrium columns are a
/// reaction-model diagnostic even when the supplied composition is frozen.
pub fn transport(state: &State) -> Result<Transport> {
    let t = state.temperature_k;
    let db = database()?;
    let gas: Vec<_> = state
        .amounts
        .iter()
        .copied()
        .filter(|(_, s)| s.data().phase() == 0)
        .collect();
    let ng = gas.len();
    let total_gas: f64 = gas.iter().map(|(n, _)| n).sum();
    if ng == 0 || total_gas <= 0.0 {
        return Ok(Transport {
            cp_frozen_j_kg_k: state.cp_frozen_j_kg_k,
            cp_equilibrium_j_kg_k: state.cp_equilibrium_j_kg_k,
            ..Transport::default()
        });
    }
    let (components, a) = basis(&gas, state.transport_atomic)?;
    let mut selected = components.clone();
    let mut total: f64 = selected.iter().map(|&i| gas[i].0).sum();
    let mut threshold = total_gas / ng as f64;
    for _ in 0..ng {
        if total >= 0.999999999 * total_gas * (1.0 - 1e-6) || selected.len() >= 40 {
            break;
        }
        threshold /= 10.0;
        for (i, &(n, _)) in gas.iter().enumerate() {
            if n >= threshold && !selected.contains(&i) && selected.len() < 40 {
                selected.push(i);
                total += n;
            }
        }
        if threshold < 1e-11 * total_gas {
            break;
        }
    }
    let nm = selected.len();
    let species: Vec<_> = selected.iter().map(|&i| gas[i].1).collect();
    let xs: Vec<_> = selected.iter().map(|&i| gas[i].0 / total).collect();
    let mw: Vec<_> = species.iter().map(|s| s.data().mw()).collect();
    let cp: Vec<_> = species.iter().map(|s| s.data().cp_over_r(t)).collect();
    let charge: Vec<_> = species
        .iter()
        .map(|&s| coefficient(s, Constituent::E))
        .collect();
    let ions = gas
        .iter()
        .any(|(_, s)| coefficient(*s, Constituent::E) != 0.0);
    let mut eta = DMatrix::zeros(nm, nm);
    let mut cond = vec![0.0; nm];
    let pure: Vec<_> = species
        .iter()
        .map(|s| {
            db.iter()
                .find(|r| r.names[0] == s.data().symbol() && r.names[1].is_empty())
        })
        .collect();
    for i in 0..nm {
        if let Some(record) = pure[i] {
            eta[(i, i)] = evaluate(&record.viscosity, t);
            cond[i] = evaluate(&record.conductivity, t);
        }
        for j in 0..i {
            if pure[i].is_some()
                && pure[j].is_some()
                && let Some(record) = db.iter().find(|r| {
                    (r.names[0] == species[i].data().symbol()
                        && r.names[1] == species[j].data().symbol())
                        || (r.names[0] == species[j].data().symbol()
                            && r.names[1] == species[i].data().symbol())
                })
            {
                eta[(i, j)] = evaluate(&record.viscosity, t);
                eta[(j, i)] = eta[(i, j)];
            }
        }
    }
    const BOLTZ: f64 = 1.380658;
    const AVOGADRO: f64 = 6.0221367;
    let viscns = 0.3125 * (1e5 * BOLTZ / (PI * AVOGADRO)).sqrt();
    let ekt = 4.8032_f64.powi(2) / (BOLTZ * t / 1000.0);
    let qc = 100.0 * ekt.powi(2);
    let xel = species
        .iter()
        .position(|s| s.data().mw() < 1.0)
        .map_or(0.0, |i| xs[i])
        .max(1e-12);
    let debye = (22.5 / PI) * (R / AVOGADRO * 100.0) * (t / 1000.0 / xel) / ekt.powi(3);
    let ionic = ((810.0 / (4.0 * PI)) * (R / AVOGADRO * 100.0) * (t / 1000.0 / xel))
        .powf(2.0 / 3.0)
        / ekt.powi(2);
    #[allow(clippy::approx_constant)] // Preserve NASA's rounded Coulomb-log floor.
    let lambda = (debye + ionic).sqrt().max(2.71828183);
    for i in 0..nm {
        if ions && (charge[i].abs() - 1.0).abs() < 1e-8 && eta[(i, i)].abs() < 1e-8 {
            continue;
        }
        if eta[(i, i)].abs() < 1e-8 {
            let omega = (50.0 * mw[i].powf(4.6) / t.powf(1.4)).ln().max(1.0);
            eta[(i, i)] = viscns * (mw[i] * t).sqrt() / omega;
        }
        if cond[i].abs() < 1e-8 {
            cond[i] = eta[(i, i)] * R * (0.00375 + 0.00132 * (cp[i] - 2.5)) / mw[i];
        }
    }
    for i in 0..nm {
        for j in i..nm {
            if eta[(i, j)].abs() >= 1e-8 {
                continue;
            }
            let mut omega = 0.0;
            if ions {
                let (ion1, ion2) = (charge[i].abs() == 1.0, charge[j].abs() == 1.0);
                if ion1 && ion2 {
                    omega = 1.36 * qc * lambda.ln();
                }
                if (ion1 && mw[j] < 1.0) || (ion2 && mw[i] < 1.0) {
                    omega = 1.29 * qc * lambda.ln();
                }
                if ion1 != ion2 {
                    omega = (6.776 - 0.4 * t.ln()).exp();
                }
            }
            if omega.abs() > 1e-8 {
                eta[(i, j)] =
                    viscns * (2.0 * t * mw[i] * mw[j] / (mw[i] + mw[j])).sqrt() * PI / omega;
                if i == j {
                    cond[i] = eta[(i, i)] * R * (0.00375 + 0.00132 * (cp[i] - 2.5)) / mw[i];
                }
            } else {
                let ratio = (mw[j] / mw[i]).sqrt();
                eta[(i, j)] = 5.656854 * eta[(i, i)] * (mw[j] / (mw[i] + mw[j])).sqrt()
                    / (1.0 + (ratio * eta[(i, i)] / eta[(j, j)]).sqrt()).powi(2);
            }
            eta[(j, i)] = eta[(i, j)];
        }
    }
    let mut viscosity = 0.0;
    let mut frozen = 0.0;
    for i in 0..nm {
        let (mut sumv, mut sumc) = (xs[i], xs[i]);
        for j in 0..nm {
            if j == i {
                continue;
            }
            let phi = 2.0 * mw[j] * eta[(i, i)] / (eta[(i, j)] * (mw[i] + mw[j]));
            let psi = phi
                * (1.0
                    + 2.41 * (mw[i] - mw[j]) * (mw[i] - 0.142 * mw[j]) / (mw[i] + mw[j]).powi(2));
            sumv += phi * xs[j];
            sumc += psi * xs[j];
        }
        viscosity += eta[(i, i)] * xs[i] / sumv;
        frozen += cond[i] * xs[i] / sumc;
    }
    let (cp_reaction, conductivity_reaction) =
        reactions(&species, &selected, &components, &a, &xs, &mw, &eta, t)?;
    let weight: f64 = xs.iter().zip(&mw).map(|(x, m)| x * m).sum();
    let cp_frozen = R * xs.iter().zip(&cp).map(|(x, c)| x * c).sum::<f64>() / weight;
    let cp_equilibrium = cp_frozen + cp_reaction / weight * 1000.0;
    // NASA stores millipoise and mW/(cm K); convert both to SI here.
    viscosity *= 1e-7;
    frozen *= 1e-4;
    let equilibrium = frozen + conductivity_reaction * 1e-4;
    let result = Transport {
        viscosity_pa_s: viscosity,
        cp_frozen_j_kg_k: cp_frozen,
        cp_equilibrium_j_kg_k: cp_equilibrium,
        conductivity_frozen_w_m_k: frozen,
        conductivity_equilibrium_w_m_k: equilibrium,
        prandtl_frozen: viscosity * cp_frozen / frozen,
        prandtl_equilibrium: viscosity * cp_equilibrium / equilibrium,
    };
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn reactions(
    species: &[Species],
    selected: &[usize],
    components: &[usize],
    a: &DMatrix<f64>,
    xs: &[f64],
    mw: &[f64],
    eta: &DMatrix<f64>,
    t: f64,
) -> Result<(f64, f64)> {
    let nm = species.len();
    let mut local: Vec<_> = components
        .iter()
        .filter_map(|j| selected.iter().position(|i| i == j))
        .collect();
    if let Some(electron) = species.iter().position(|s| s.data().mw() < 1.0)
        && !local.contains(&electron)
    {
        local.push(electron);
    }
    let mut alpha = Vec::new();
    for i in 0..nm {
        if local.contains(&i) {
            continue;
        }
        let mut row = DVector::zeros(nm);
        row[i] = -1.0;
        for (r, &j) in local.iter().enumerate() {
            row[j] = if r < components.len() {
                a[(r, selected[i])]
            } else {
                coefficient(species[i], Constituent::E)
            };
        }
        alpha.push(row);
    }
    // Eliminate negligible species from all but one reaction, as NASA does.
    for i in 0..nm {
        if xs[i] >= 1e-10 {
            continue;
        }
        let mut first: Option<DVector<f64>> = None;
        for row in &mut alpha {
            let coeff = row[i];
            if coeff.abs() > 1e-5 {
                if let Some(first) = &first {
                    *row = &*row / coeff - first;
                } else {
                    first = Some(&*row / coeff);
                }
            }
        }
    }
    let nr = alpha.len();
    if nr == 0 {
        return Ok((0.0, 0.0));
    }
    let dh = DVector::from_iterator(
        nr,
        alpha.iter().map(|row| {
            (0..nm)
                .map(|j| row[j] * species[j].data().h_over_r(t) / t)
                .sum()
        }),
    );
    for row in &mut alpha {
        for x in row.iter_mut() {
            if x.abs() < 1e-6 {
                *x = 0.0;
            }
        }
    }
    let mut g = DMatrix::zeros(nr, nr);
    let mut diffusion = DMatrix::zeros(nr, nr);
    for k in 0..nm {
        for m in k + 1..nm {
            if xs[k] < 1e-10 || xs[m] < 1e-10 {
                continue;
            }
            let resistance = mw[k] * mw[m] / (1.1 * eta[(k, m)] * (mw[k] + mw[m]));
            let st =
                DVector::from_iterator(nr, alpha.iter().map(|row| xs[m] * row[k] - xs[k] * row[m]));
            let contribution = &st * st.transpose() / (xs[k] * xs[m]);
            diffusion += &contribution * resistance;
            g += contribution;
        }
    }
    let solve = |matrix: DMatrix<f64>| {
        let mut augmented = DMatrix::zeros(nr, nr + 1);
        augmented.view_mut((0, 0), (nr, nr)).copy_from(&matrix);
        augmented.column_mut(nr).copy_from(&dh);
        // NASA uses the overwritten right-hand side even when elimination fails.
        let _ = crate::equilibrium_core::gauss(&mut augmented);
        (R * 1e-3) * dh.dot(&augmented.column(nr))
    };
    Ok((solve(g), 0.6 * solve(diffusion)))
}

impl State {
    pub fn transport(&self) -> Result<Transport> {
        transport(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malformed_transport_records_are_data_errors() {
        let text = include_str!("../vendor/nasa-transport/trans.inp");
        for malformed in [
            text.replacen("V3C3", "V4C3", 1),
            text.replacen("0.61205763E 00", "          NaN", 1),
            text.replacen("1000.0", " 100.0", 1),
            text.replacen(" V  200.0", " X  200.0", 1),
            text[..text.len() / 2].into(),
            text.replace("end", ""),
            format!("{text}trailing"),
        ] {
            assert!(matches!(parse(&malformed), Err(Error::Data(_))));
        }
    }
}
