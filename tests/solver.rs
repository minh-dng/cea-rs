use cea_rs::{
    Error,
    rocket::FlowMode,
    solver::{Mode, Solver},
    thermo::Reactant,
};
fn hydrogen() -> Solver {
    Solver::new(
        vec![
            Reactant::named("H2").unwrap(),
            Reactant::named("O2").unwrap(),
        ],
        false,
        false,
        None,
        vec![],
    )
    .unwrap()
}
#[test]
fn six_modes_conserve_same_state() {
    let s = hydrogen();
    let w = [1.0, 8.0];
    let tp = s.solve(Mode::TP(3000.0, 1e6), &w).unwrap();
    let modes = [
        Mode::HP(tp.enthalpy_j_kg, 1e6),
        Mode::SP(tp.entropy_j_kg_k, 1e6),
        Mode::TV(3000.0, 1.0 / tp.density_kg_m3),
        Mode::UV(tp.energy_j_kg, 1.0 / tp.density_kg_m3),
        Mode::SV(tp.entropy_j_kg_k, 1.0 / tp.density_kg_m3),
    ];
    for m in modes {
        let st = s.solve(m, &w).unwrap();
        assert!((st.temperature_k / 3000.0 - 1.0).abs() < 3e-5);
        assert!(st.element_residual < 2e-5);
    }
}
#[test]
fn rejects_bad_weights_and_domain() {
    let s = hydrogen();
    for w in [
        vec![],
        vec![0.0, 0.0],
        vec![-1.0, 1.0],
        vec![f64::NAN, 1.0],
        vec![f64::INFINITY, 1.0],
    ] {
        assert!(s.solve(Mode::TP(3000.0, 1e5), &w).is_err());
    }
    for mode in [
        Mode::TP(0.0, 1e5),
        Mode::TP(f64::NAN, 1e5),
        Mode::TP(3000.0, -1.0),
        Mode::TV(3000.0, 0.0),
    ] {
        assert!(s.solve(mode, &[1.0, 8.0]).is_err());
    }
    assert!(matches!(
        Reactant::named("RP-1").unwrap().properties(500.0),
        Err(Error::Domain(_))
    ));
}
#[test]
fn graphite_forms_and_charge_balances() {
    let s = Solver::new(
        vec![Reactant::named("CO").unwrap()],
        true,
        false,
        None,
        vec![],
    )
    .unwrap();
    let st = s.solve(Mode::TP(800.0, 1e5), &[1.0]).unwrap();
    assert!(st.mole_fractions["C(gr)"] > 0.01);
    let s = Solver::new(
        vec![Reactant::named("H2").unwrap()],
        false,
        true,
        Some(vec!["H".into(), "H2".into(), "H+".into(), "e-".into()]),
        vec![],
    )
    .unwrap();
    let st = s.solve(Mode::TP(10000.0, 1e5), &[1.0]).unwrap();
    assert!(st.mole_fractions["e-"] > 1e-4);
    assert!(st.element_residual < 2e-5);
}
#[test]
fn rocket_closes_mass_entropy_energy_for_three_flow_modes() {
    let s = hydrogen();
    for flow in [
        FlowMode::Equilibrium,
        FlowMode::FrozenAtChamber,
        FlowMode::FrozenAtThroat,
    ] {
        let r = s
            .rocket(&[1.0, 6.0], &[298.15, 298.15], 3e6, 10.0, flow)
            .unwrap();
        assert!((r.throat.mach - 1.0).abs() < 1e-4);
        assert!(r.exit.mach > 1.0);
        let gt = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
        let ge = r.exit.state.density_kg_m3 * r.exit.velocity_m_s * 10.0;
        assert!((gt / ge - 1.0).abs() < 1e-4);
        assert!((r.exit.state.entropy_j_kg_k / r.chamber.state.entropy_j_kg_k - 1.0).abs() < 2e-5);
        assert!(r.cstar_m_s > 1000.0 && r.vacuum_velocity_m_s > 3000.0);
    }
}
#[test]
fn impossible_product_pool_fails_without_panicking() {
    let s = Solver::new(
        vec![
            Reactant::named("H2").unwrap(),
            Reactant::named("O2").unwrap(),
        ],
        false,
        false,
        Some(vec!["H2O".into()]),
        vec![],
    )
    .unwrap();
    assert!(matches!(
        s.solve(Mode::TP(3000.0, 1e5), &[1.0, 1.0]),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn nozzle_search_keeps_closure_near_throat_and_at_large_area() {
    let s = hydrogen();
    for flow in [
        FlowMode::Equilibrium,
        FlowMode::FrozenAtChamber,
        FlowMode::FrozenAtThroat,
    ] {
        for area in [1.001, 100.0] {
            let r = s
                .rocket(&[1.0, 6.0], &[298.15, 298.15], 3e6, area, flow)
                .unwrap();
            let throat_flux = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
            let exit_flux = r.exit.state.density_kg_m3 * r.exit.velocity_m_s * area;
            assert!((throat_flux / exit_flux - 1.0).abs() < 1e-4);
            assert!(r.exit.state.pressure_pa <= r.throat.state.pressure_pa * (1.0 + 1e-8));
            assert!(r.exit.velocity_m_s.is_finite());
        }
    }
}

#[test]
fn carbon_nozzle_preserves_chamber_phase_initialization() {
    let reactant = |name: &str, formula: &[(&str, f64)], mw, enthalpy| {
        Reactant::custom(
            name.into(),
            formula.iter().map(|(e, n)| ((*e).into(), *n)).collect(),
            mw,
            enthalpy,
            298.15,
        )
        .unwrap()
    };
    let s = Solver::new(
        vec![
            reactant(
                "C3H3N_reference",
                &[("C", 3.0), ("H", 3.0), ("N", 1.0)],
                53.0626,
                140000.0,
            ),
            reactant(
                "C4H6_reference",
                &[("C", 4.0), ("H", 6.0)],
                54.0904,
                90500.0,
            ),
            reactant(
                "C8H8_reference",
                &[("C", 8.0), ("H", 8.0)],
                104.1491,
                103400.0,
            ),
            reactant("N2O_ref298", &[("N", 2.0), ("O", 1.0)], 44.0128, 81600.0),
        ],
        true,
        false,
        None,
        vec![],
    )
    .unwrap();
    // Regression point from a carbon-bearing Hybrid alpha.2 burn.
    let weights = [
        0.2015848925668829,
        0.23920531966165437,
        0.06418196354020744,
        0.49502782423125524,
    ];
    let pressure = 417101.1025245694;
    for flow in [FlowMode::Equilibrium, FlowMode::FrozenAtThroat] {
        for area in [7.5608953312955025] {
            let r = s
                .rocket(&weights, &[298.15; 4], pressure, area, flow)
                .unwrap();
            assert!(r.chamber.state.mole_fractions["C(gr)"] > 0.0);
            assert!((r.throat.mach - 1.0).abs() < 1e-4);
            let flux = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
            assert!(
                (flux / (r.exit.state.density_kg_m3 * r.exit.velocity_m_s * area) - 1.0).abs()
                    < 1e-4
            );
            for st in [&r.chamber, &r.throat, &r.exit] {
                assert!(st.state.element_residual < 2e-5);
                assert!(st.state.constraint_residual < 2e-5);
                assert!(
                    (st.state.entropy_j_kg_k / r.chamber.state.entropy_j_kg_k - 1.0).abs() < 2e-5
                );
                assert!(
                    ((st.state.enthalpy_j_kg + st.velocity_m_s.powi(2) / 2.0)
                        / r.chamber.state.enthalpy_j_kg
                        - 1.0)
                        .abs()
                        < 2e-5
                );
            }
            assert!(r.cstar_m_s.is_finite() && r.cstar_m_s > 0.0);
            assert!(r.vacuum_velocity_m_s.is_finite() && r.vacuum_velocity_m_s > 0.0);
            assert!(r.exit.mach > 1.0);
        }
    }
}

fn fixture_reactant(name: &str) -> Reactant {
    if name == "CHOS-Binder" {
        Reactant::custom(
            name.into(),
            [("C", 1.0), ("H", 1.86955), ("O", 0.031256), ("S", 0.008415)]
                .into_iter()
                .map(|(e, n)| (e.into(), n))
                .collect(),
            14.6652984484,
            -2999.082 * 4.184,
            298.15,
        )
        .unwrap()
    } else {
        Reactant::named(name).unwrap()
    }
}

#[test]
fn native_equilibrium_domain_regressions() {
    // Oracle and regeneration command are in scripts/record_equilibrium_parity.py.
    let numbers = |s: &str| {
        s.split(';')
            .filter(|x| !x.is_empty())
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<_>>()
    };
    let names = |s: &str| {
        s.split(';')
            .filter(|x| !x.is_empty())
            .map(String::from)
            .collect::<Vec<_>>()
    };
    let mut worst = 0.0_f64;
    let mut worst_x = 0.0_f64;
    for line in include_str!("fixtures/equilibrium-native.txt")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let c: Vec<_> = line.split('|').collect();
        let weights = numbers(c[2]);
        let solver = Solver::new(
            c[1].split(';').map(fixture_reactant).collect(),
            true,
            false,
            if c[7].is_empty() {
                None
            } else {
                Some(names(c[7]))
            },
            vec![],
        )
        .unwrap()
        .with_omit(names(c[8]))
        .unwrap()
        .with_trace(c[6].parse().unwrap())
        .unwrap();
        let first = c[4].parse().unwrap();
        let second = c[5].parse().unwrap();
        let mode = match c[3] {
            "TP" => Mode::TP(first, second),
            "HP" => Mode::HP(first, second),
            "SP" => Mode::SP(first, second),
            "TV" => Mode::TV(first, second),
            "UV" => Mode::UV(first, second),
            _ => panic!("Unknown fixture mode"),
        };
        if !c[9].is_empty() {
            let h = solver.reactant_enthalpy(&weights, &numbers(c[9])).unwrap();
            assert!(
                (h - first).abs() <= 1e-8 * first.abs().max(1.0),
                "{}: reactant enthalpy {h} != {first}",
                c[0]
            );
        }
        let state = solver
            .solve(mode, &weights)
            .unwrap_or_else(|e| panic!("{}: {e}", c[0]));
        let values = [
            state.temperature_k,
            state.pressure_pa,
            state.density_kg_m3,
            state.enthalpy_j_kg,
            state.entropy_j_kg_k,
            state.cp_equilibrium_j_kg_k,
            state.gamma_s,
            state.molecular_weight_kg_kmol,
        ];
        for (i, (actual, expected)) in values.into_iter().zip(numbers(c[10])).enumerate() {
            let difference = if expected == 0.0 {
                actual.abs()
            } else {
                ((actual - expected) / expected).abs()
            };
            worst = worst.max(difference);
            assert!(
                difference <= 1e-4,
                "{} field {i}: {actual} != {expected} (relative {difference})",
                c[0]
            );
        }
        let expected: std::collections::BTreeMap<_, _> = c[11]
            .split(';')
            .map(|pair| {
                let (name, value) = pair.split_once('=').unwrap();
                (name.to_string(), value.parse::<f64>().unwrap())
            })
            .collect();
        for name in expected.keys().chain(state.mole_fractions.keys()) {
            let difference = (expected.get(name).unwrap_or(&0.0)
                - state.mole_fractions.get(name).unwrap_or(&0.0))
            .abs();
            worst_x = worst_x.max(difference);
            assert!(
                difference <= 3e-5,
                "{} mole fraction {name}: difference {difference}",
                c[0]
            );
        }
    }
    println!(
        "Native equilibrium goldens: worst relative property difference={worst:e}, absolute mole fraction difference={worst_x:e}"
    );
}

#[test]
fn native_reactant_reference_enthalpies_and_restored_records() {
    // cea==3.3.4 Mixture([name]).calc_property(ENTHALPY, [1.0], T), in J/kg.
    for (name, t, expected) in [
        ("N2O", 298.15, 1854006.1071421793),
        ("N2O", 250.0, 1812994.569278336),
        ("Air", 298.15, -4333.823194800017),
        ("Air", 7000.0, 0.0),
        ("AL(cr)", 298.15, 0.0015293115434638815),
        ("CH4(L)", 111.64, -5562301.54228217),
        ("RP-1", 500.0, -1768558.6973210066),
    ] {
        let solver = Solver::new(
            vec![Reactant::named(name).unwrap()],
            true,
            false,
            None,
            vec![],
        )
        .unwrap();
        let h = solver.reactant_enthalpy(&[1.0], &[t]).unwrap();
        assert!(
            (h - expected).abs() <= 1e-7 * expected.abs().max(1.0),
            "{name}: {h} != {expected}"
        );
    }
    let methane = Reactant::named("CH4(L)").unwrap();
    assert!(methane.properties(111.64).is_ok());
    assert!(
        methane
            .properties(methane.reference_temperature + 10.0)
            .is_ok()
    );
    assert!(
        methane
            .properties(methane.reference_temperature + 10.01)
            .is_err()
    );

    // Same native call at 350 K, converted to J/mol by multiplying M/1000.
    for (name, expected) in [
        ("Br2(cr)", -7204.888557889754),
        ("Ca(a)", 1351.2922740037761),
        ("CrN(cr)", -114675.50761929515),
        ("FeCL3(cr)", -394097.66498064925),
        ("FeOCL(cr)", -407256.6824007526),
        ("Fe3O4(cr)", -1110183.0954285394),
        ("Li(cr)", 1312.6528028556324),
        ("NH4F(cr)", -463940.1652045158),
        ("Si(cr)", 1068.0958048973955),
        ("Ti3O5(a)", -2447778.298560001),
        ("U3O8(II)", -3561923.662937205),
    ] {
        let h = Reactant::named(name).unwrap().properties(350.0).unwrap().1;
        assert!(
            (h / expected - 1.0).abs() < 1e-8,
            "{name}: {h} != {expected}"
        );
    }
}

#[test]
fn omission_trace_and_extrapolation_are_explicit() {
    let solver = hydrogen()
        .with_omit(vec!["H2O".into()])
        .unwrap()
        .with_trace(1e-15)
        .unwrap();
    let state = solver.solve(Mode::TP(3000.0, 1e6), &[1.0, 8.0]).unwrap();
    assert!(!state.mole_fractions.contains_key("H2O"));
    assert!(state.element_residual < 2e-5);
    let low = hydrogen().solve(Mode::TP(180.0, 1e5), &[1.0, 8.0]).unwrap();
    assert!(low.extrapolated_species.contains(&"H2O".to_string()));
    assert!(
        low.trace_extrapolated_species
            .iter()
            .all(|n| low.mole_fractions[n] <= 1e-11)
    );
    let high = hydrogen()
        .solve(Mode::TP(10000.0, 1e6), &[1.0, 8.0])
        .unwrap();
    assert!(high.extrapolated_species.contains(&"H2O".to_string()));
    assert!(hydrogen().with_trace(f64::NAN).is_err());
    assert!(
        hydrogen()
            .with_omit(vec!["unknown-product".into()])
            .unwrap()
            .solve(Mode::TP(3000.0, 1e6), &[1.0, 8.0])
            .unwrap()
            .mole_fractions
            .contains_key("H2O")
    );
}
