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
        assert!((r.throat.mach - 1.0).abs() < 2e-6);
        assert!(r.exit.mach > 1.0);
        let gt = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
        let ge = r.exit.state.density_kg_m3 * r.exit.velocity_m_s * 10.0;
        assert!((gt / ge - 1.0).abs() < 2e-6);
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
    assert!(s.solve(Mode::TP(3000.0, 1e5), &[1.0, 1.0]).is_err());
}

#[test]
fn nozzle_search_keeps_closure_near_throat_and_at_large_area() {
    let s = hydrogen();
    for flow in [
        FlowMode::Equilibrium,
        FlowMode::FrozenAtChamber,
        FlowMode::FrozenAtThroat,
    ] {
        for area in [1.0, 1.001, 100.0] {
            let r = s
                .rocket(&[1.0, 6.0], &[298.15, 298.15], 3e6, area, flow)
                .unwrap();
            let throat_flux = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
            let exit_flux = r.exit.state.density_kg_m3 * r.exit.velocity_m_s * area;
            assert!((throat_flux / exit_flux - 1.0).abs() < 2e-6);
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
    let weights = [
        0.2015848925668829,
        0.23920531966165437,
        0.06418196354020744,
        0.49502782423125524,
    ];
    let pressure = 417101.1025245694;
    for area in [1.0, 7.5608953312955025] {
        let r = s
            .rocket(
                &weights,
                &[298.15; 4],
                pressure,
                area,
                FlowMode::Equilibrium,
            )
            .unwrap();
        assert!(r.chamber.state.mole_fractions["C(gr)"] > 0.0);
        assert!((r.throat.mach - 1.0).abs() < 2e-6);
        let flux = r.throat.state.density_kg_m3 * r.throat.velocity_m_s;
        assert!(
            (flux / (r.exit.state.density_kg_m3 * r.exit.velocity_m_s * area) - 1.0).abs() < 2e-6
        );
        for st in [&r.chamber, &r.throat, &r.exit] {
            assert!(st.state.equilibrium);
            assert!(st.state.element_residual < 2e-5);
            assert!(st.state.constraint_residual < 2e-5);
            assert!((st.state.entropy_j_kg_k / r.chamber.state.entropy_j_kg_k - 1.0).abs() < 2e-5);
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
        if area > 1.0 {
            assert!(r.exit.mach > 1.0);
        }
    }
}
