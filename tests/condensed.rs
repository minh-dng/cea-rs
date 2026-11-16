use cea_rs::{
    R,
    solver::{Mode, Solver, State},
    thermo::Reactant,
};
use thermo_species::Species;

fn model(name: &str) -> Solver {
    Solver::new(
        vec![
            Reactant::named("Ar").unwrap(),
            Reactant::named(name).unwrap(),
        ],
        true,
        false,
        None,
        vec![],
    )
    .unwrap()
}

// Recompute balances and chemical potentials from the reported composition.
fn check(
    state: &State,
    solver: &Solver,
    weights: &[f64],
    vapor: &str,
    phases: &[&str],
) -> (f64, f64) {
    let amounts: Vec<_> = state
        .mass_fractions
        .iter()
        .filter(|(_, y)| **y > 0.0)
        .map(|(name, y)| {
            let sp = name.parse::<Species>().unwrap();
            (y / sp.data().mw(), sp)
        })
        .collect();
    let ng: f64 = amounts
        .iter()
        .filter(|(_, sp)| sp.data().phase() == 0)
        .map(|(n, _)| n)
        .sum();
    let t = state.temperature_k;
    let mut h = 0.0;
    let mut s = 0.0;
    for &(n, sp) in &amounts {
        let d = sp.data();
        h += n * d.h(t) * 1000.0;
        s += n
            * (d.s(t) * 1000.0
                - if d.phase() == 0 {
                    R * (n / ng * state.pressure_pa / 1e5).ln()
                } else {
                    0.0
                });
    }
    for element in solver
        .reactants
        .iter()
        .flat_map(|r| &r.formula)
        .map(|(_, e)| e)
    {
        let feed: f64 = solver
            .reactants
            .iter()
            .zip(weights)
            .map(|(r, w)| {
                w / weights.iter().sum::<f64>() / r.molecular_weight
                    * r.formula
                        .iter()
                        .filter(|(_, e)| e == element)
                        .map(|(a, _)| a)
                        .sum::<f64>()
            })
            .sum();
        let product: f64 = amounts
            .iter()
            .map(|(n, sp)| {
                n * sp
                    .data()
                    .constituents()
                    .iter()
                    .filter(|(_, e)| e == element)
                    .map(|(a, _)| a)
                    .sum::<f64>()
            })
            .sum();
        assert!(
            (product / feed - 1.0).abs() < 2e-6,
            "{element:?}: {product} != {feed}"
        );
    }
    assert!((h - state.enthalpy_j_kg).abs() < 0.02, "h: {h}");
    assert!((s - state.entropy_j_kg_k).abs() < 1e-4, "s: {s}");
    let mu = |name: &str| {
        let sp = name.parse::<Species>().unwrap();
        (sp.data().h(t) / t - sp.data().s(t)) * 1000.0 / R
    };
    let lowest = phases
        .iter()
        .map(|name| mu(name))
        .fold(f64::INFINITY, f64::min);
    let mut present = 0;
    for &name in phases {
        let gap = mu(name) - lowest;
        if *state.mass_fractions.get(name).unwrap_or(&0.0) > 0.0 {
            present += 1;
            assert!(gap.abs() < 1e-7, "{name}: mu gap={gap}");
        }
    }
    assert!(present > 0);
    let nv = amounts
        .iter()
        .filter(|(_, sp)| sp.data().symbol() == vapor)
        .map(|(n, _)| n)
        .sum::<f64>();
    let supersaturation = mu(vapor) + (nv / ng * state.pressure_pa / 1e5).ln() - lowest;
    assert!(
        supersaturation <= 1e-7,
        "{vapor}: ln(p/p_sat)={supersaturation}"
    );
    if nv / ng > 1e-11 {
        assert!(supersaturation.abs() < 1e-7);
    }
    (h, s)
}

#[test]
fn unseeded_cold_enthalpy_and_entropy_find_condensed_equilibrium() {
    let solver = model("Hg(L)");
    let weights = [1.0, 20.0];
    let target_h = -15551.803948663037;
    let a = solver.solve(Mode::TP(232.0, 1e5), &weights).unwrap();
    let b = solver.solve(Mode::TP(237.0, 1e5), &weights).unwrap();
    let target_s = (a.entropy_j_kg_k + b.entropy_j_kg_k) / 2.0;
    for mode in [Mode::HP(target_h, 1e5), Mode::SP(target_s, 1e5)] {
        let state = solver.solve(mode, &weights).unwrap();
        let (h, s) = check(&state, &solver, &weights, "Hg", &["Hg(cr)", "Hg(L)"]);
        assert!(state.mass_fractions["Hg(cr)"] > 0.1);
        assert!(state.mass_fractions["Hg(L)"] > 0.1);
        match mode {
            Mode::HP(..) => assert!((h - target_h).abs() < 0.02),
            _ => assert!((s - target_s).abs() < 1e-4),
        }
    }
}

#[test]
fn reversed_condensed_interval_does_not_hide_supersaturated_bromine() {
    let solver = model("Br2(L)");
    let state = solver.solve(Mode::TP(262.0, 1e5), &[1.0, 20.0]).unwrap();
    assert_eq!(state.temperature_k, 262.0);
    check(&state, &solver, &[1.0, 20.0], "Br2", &["Br2(cr)", "Br2(L)"]);
    assert!(state.mass_fractions["Br2(cr)"] > 0.1);
}
