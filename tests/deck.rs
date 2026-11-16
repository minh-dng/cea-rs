use cea_rs::{
    R,
    deck::{AmountBasis, Problem, Ratio, parse_deck},
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
",
    )
    .unwrap();
    assert_eq!(decks.len(), 4);
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
