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
            // At a tabulated crossing both phases are equal within the
            // rounding of that crossing temperature (about 1e-8 K).
            assert!(gap.abs() < 1e-9, "{name}: mu gap={gap}");
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

#[test]
fn fixed_temperature_selects_the_lower_gibbs_phase_without_coexistence() {
    for (reactant, t, vapor, phases) in [
        ("H2O(L)", 273.15, "H2O", ["H2O(cr)", "H2O(L)"]),
        ("H2O(L)", 273.121944434780, "H2O", ["H2O(cr)", "H2O(L)"]),
        ("Hg(L)", 234.29, "Hg", ["Hg(cr)", "Hg(L)"]),
        ("Hg(L)", 234.28999933688027, "Hg", ["Hg(cr)", "Hg(L)"]),
        ("NaCN", 288.5, "NaCN", ["NaCN(II)", "NaCN(III)"]),
        ("NaCN", 835.0, "NaCN", ["NaCN(III)", "NaCN(L)"]),
    ] {
        let solver = model(reactant);
        let tp = solver.solve(Mode::TP(t, 1e5), &[1.0, 1.0]).unwrap();
        let tv = solver
            .solve(Mode::TV(t, 1.0 / tp.density_kg_m3), &[1.0, 1.0])
            .unwrap();
        for state in [tp, tv] {
            assert_eq!(state.temperature_k, t);
            check(&state, &solver, &[1.0, 1.0], vapor, &phases);
            assert_eq!(
                phases
                    .iter()
                    .filter(|&&name| *state.mass_fractions.get(name).unwrap_or(&0.0) > 0.0)
                    .count(),
                1
            );
        }
    }
    // Below the lambda transition, NaCN(III) is used outside its fit and must not
    // replace NaCN(II), although its extrapolated Gibbs energy is lower.
    let solver = model("NaCN");
    for t in [250.0, 285.0] {
        let state = solver.solve(Mode::TP(t, 1e5), &[1.0, 1.0]).unwrap();
        assert!(state.mass_fractions["NaCN(II)"] > 0.4);
        assert!(
            !state.mass_fractions.contains_key("NaCN(III)")
                || state.mass_fractions["NaCN(III)"] == 0.0
        );
    }
}

#[test]
fn nacn_at_300_k_selects_the_phase_inside_its_fit() {
    for reactant in ["NaCN(II)", "NaCN"] {
        let solver = model(reactant);
        let state = solver.solve(Mode::TP(300.0, 1e5), &[1.0, 1.0]).unwrap();
        assert_eq!(state.temperature_k, 300.0);
        check(&state, &solver, &[1.0, 1.0], "NaCN", &["NaCN(III)"]);
        assert!(state.mass_fractions["NaCN(III)"] > 0.49);
        assert!(!state.mass_fractions.contains_key("NaCN(II)"));
        assert!(!state.mass_fractions.contains_key("NaCN(L)"));
    }
}

#[test]
fn seeded_nacn_changes_phase_at_fixed_temperature() {
    let solver = Solver::new(
        model("NaCN").reactants,
        true,
        false,
        None,
        vec!["NaCN(II)".into()],
    )
    .unwrap();
    for (t, phases, expected) in [
        (300.0, vec!["NaCN(III)"], "NaCN(III)"),
        (835.0, vec!["NaCN(III)", "NaCN(L)"], "NaCN(L)"),
    ] {
        let state = solver.solve(Mode::TP(t, 1e5), &[1.0, 1.0]).unwrap();
        assert_eq!(state.temperature_k, t);
        check(&state, &solver, &[1.0, 1.0], "NaCN", &phases);
        assert!(state.mass_fractions[expected] > 0.49);
        let (low, high) = expected
            .parse::<Species>()
            .unwrap()
            .data()
            .valid_temperature_range();
        assert!((low..=high).contains(&state.temperature_k));
        assert!(!state.mass_fractions.contains_key("NaCN(II)"));
    }
}
