use cea_rs::{
    R,
    deck::{AmountBasis, CaseInput, CaseOutput, Problem, Ratio, parse_deck, run_deck},
    solver::{Mode, Solver},
    thermo::Reactant,
};

#[test]
fn parses_units_formulas_punctuation_and_multiple_decks() {
    let decks = parse_deck(
        "\
reac name Binder C 1 H 1.86955 O .031256 S .008415 wt%=100 h,cal=-2999.082 t,c=25
problem hp case=custom p,psia=500,250
output calories transport long
end trailing text is a NASA end-line comment
problem uv u/r=-45.1343 rho,kg/m**3=14.428
reac name C8H18(L),n-octa moles=1
end
problem tp t,f=32,212 p,mmhg=10
reac name H2 moles=1
end
problem tp t,r=540 p,atm=1
reac name H2
end
problem hp p=1
reac name Air h,j=0 t=298.15
end
problem hp p=1
reac name Heavy H 2 mw=3 h,j=0 t=298.15
end
",
    )
    .unwrap();
    assert_eq!(decks.len(), 6);
    let d = &decks[0];
    assert_eq!(d.problem, Problem::HP);
    assert_eq!(d.pressures_pa, [500.0 * 6894.75729, 250.0 * 6894.75729]);
    let r = &d.reactants[0];
    let mw = 12.0107 + 1.86955 * 1.00794 + 0.031256 * 15.9994 + 0.008415 * 32.065;
    assert!((r.reactant.molecular_weight - mw).abs() < 1e-12);
    assert_eq!(r.temperature_k, Some(298.15));
    assert_eq!(r.reactant.properties(298.15).unwrap().1, -2999.082 * 4.184);
    assert!(d.output.transport);
    assert_eq!(decks[1].energies_j_kg, [-45.1343 * R]);
    assert_eq!(decks[1].volumes_m3_kg, [1.0 / 14.428]);
    assert_eq!(decks[1].reactants[0].reactant.name, "C8H18(L),n-octa");
    assert_eq!(decks[1].reactants[0].basis, AmountBasis::Moles);
    assert_eq!(decks[2].temperatures_k, [273.15, 373.15]);
    assert_eq!(decks[2].pressures_pa, [1333.0]);
    assert_eq!(decks[3].temperatures_k, [300.0]);
    assert_eq!(decks[3].pressures_pa, [101325.0]);
    assert_eq!(
        decks[4].reactants[0].reactant.molecular_weight,
        Reactant::named("Air").unwrap().molecular_weight
    );
    assert_eq!(decks[5].reactants[0].reactant.molecular_weight, 3.0);
}

#[test]
fn runs_six_modes_with_si_constraints_and_stream_enthalpy() {
    let solver = Solver::new(
        ["H2", "O2"].map(|n| Reactant::named(n).unwrap()).to_vec(),
        true,
        false,
        None,
        vec![],
    )
    .unwrap();
    let reference = solver.solve(Mode::TP(3000.0, 1e6), &[1.0, 8.0]).unwrap();
    let h = reference.enthalpy_j_kg / R;
    let u = reference.energy_j_kg / R;
    let s = reference.entropy_j_kg_k / R;
    let v = 1.0 / reference.density_kg_m3;
    for constraint in [
        "tp t=3000 p=10".to_string(),
        format!("hp h/r={h} p=10"),
        format!("sp s/r={s} p=10"),
        format!("tv t=3000 v,m**3/kg={v}"),
        format!("uv u/r={u} v,m**3/kg={v}"),
        format!("sv s/r={s} v,m**3/kg={v}"),
    ] {
        let results = run_deck(&format!(
            "problem {constraint}\nreac name H2 wt%=1 name O2 wt%=8\noutput transport\nend"
        ))
        .unwrap();
        let CaseOutput::Equilibrium(actual) = &results[0].output else {
            panic!("wrong result type")
        };
        assert!(
            (actual.temperature_k / 3000.0 - 1.0).abs() < 1e-7,
            "{constraint}"
        );
        assert!(
            (actual.pressure_pa / 1e6 - 1.0).abs() < 1e-7,
            "{constraint}"
        );
        assert_eq!(results[0].transport.len(), 1);
    }
    let deck = parse_deck("reac name Fuel C 1 H 4 wt%=1 h,j=-74873 t=298.15\nproblem hp p=10\nend")
        .unwrap();
    let cases = deck[0].cases().unwrap();
    let CaseInput::Equilibrium(Mode::HP(h, _)) = cases[0].input else {
        panic!("wrong mode")
    };
    assert!((h + 74873.0 * 1000.0 / 16.04246).abs() < 1e-8);
}

#[test]
fn mixture_ratios_shock_schedules_and_unsupported_modes_are_explicit() {
    let reference =
        run_deck("problem tp t=3000 p=10 o/f=7.93668273905193\nreac fuel H2 oxid O2\nend").unwrap();
    let of = reference[0].weights[1] / reference[0].weights[0];
    for ratio in [
        "r=1".into(),
        "phi=1".into(),
        format!("%fuel={}", 100.0 / (1.0 + of)),
        format!("f/o={}", 1.0 / of),
    ] {
        let deck = parse_deck(&format!(
            "problem tp t=3000 p=10 {ratio}\nreac fuel H2 oxid O2\nend"
        ))
        .unwrap();
        let weights = &deck[0].cases().unwrap()[0].weights;
        assert!(
            (weights[1] / weights[0] / (reference[0].weights[1] / reference[0].weights[0]) - 1.0)
                .abs()
                < 1e-10,
            "{ratio}: {weights:?}"
        );
    }
    let deck = parse_deck("problem shock p=1,2 t=300 u1=1100,1200,1300 incd refl eql froz\nreac name H2 moles=.05 name O2 moles=.05 name Ar moles=.9\nend").unwrap();
    let cases = deck[0].cases().unwrap();
    assert_eq!(cases.len(), 12);
    for (i, case) in cases.iter().enumerate() {
        let CaseInput::Shock {
            pressure_pa,
            incident_frozen,
            reflected_frozen,
            ..
        } = case.input
        else {
            panic!("wrong mode")
        };
        assert_eq!(pressure_pa, if i % 3 == 0 { 1e5 } else { 2e5 });
        assert_eq!(incident_frozen, i >= 6);
        assert_eq!(reflected_frozen, (i / 3) % 2 == 0);
        assert!(matches!(case.run().unwrap().output, CaseOutput::Shock(_)));
    }
    let error =
        run_deck("problem detonation frozen p=1 t=300\nreac fuel H2 oxid O2\nend").unwrap_err();
    assert!(matches!(error, cea_rs::Error::Unsupported(_)));
    assert!(error.to_string().contains("Detonation"));
}

#[test]
fn rejects_unknown_keywords_bad_units_and_incomplete_decks() {
    for (text, keyword) in [
        ("problem tp p,banana=1 t=300\nreac name H2\nend", "banana"),
        (
            "problem tp p=1 t=300 mystery=1\nreac name H2\nend",
            "mystery",
        ),
        (
            "problem tp p=1 t=300\nreac name H2 nonsense=2\nend",
            "nonsense",
        ),
        (
            "problem tp p=1 t=300\nreac name H2\noutput silent\nend",
            "silent",
        ),
        ("problem tp p=1 t=300\nreac name H2", "end"),
        ("problem tp p=1 t=nan\nreac name H2\nend", "finite"),
        ("problem tp p=1 t=300 h/r=1\nreac name H2\nend", "h/r"),
        ("problem tp,banana p=1 t=300\nreac name H2\nend", "banana"),
        ("problem hp p=1\nreac name Fuel C 1 H 4\nend", "requires h"),
        (
            "problem tp p=1 t=300\nreac name H2 wt%=1 name O2 moles=1\nend",
            "mix",
        ),
    ] {
        let error = parse_deck(text).unwrap_err().to_string();
        assert!(error.contains(keyword), "{keyword}: {error}");
    }
    let d = parse_deck(
        "problem tp p=1 t=300 r,eq.ratio=1,1.5\nreac fuel H2 moles=1 oxid Air moles=1\nend",
    )
    .unwrap();
    assert_eq!(d[0].ratios, Some((Ratio::Chemical, vec![1.0, 1.5])));
}
