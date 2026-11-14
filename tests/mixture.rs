//! Goldens: native PyPI cea==3.3.4, matching NASA thermo.inp, macOS arm64.
//! Construct cea.Mixture(names), then call chem_eq_ratio_to_of_ratio and
//! weight_eq_ratio_to_of_ratio with r=phi=0.75, of_ratio_to_weights with 2.5,
//! and weights_to_moles/moles_to_weights with [1,2,...]. Scalar ratios are
//! exactly representable: native Cython rounds them through float32.
use cea_rs::{
    mixture,
    solver::{Mode, Solver},
    thermo::Reactant,
};

fn solver(names: &[&str]) -> Solver {
    Solver::new(
        names.iter().map(|n| Reactant::named(n).unwrap()).collect(),
        true,
        false,
        None,
        vec![],
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
#[test]
#[rustfmt::skip]
fn nasa_mixture_conversions() {
    {
        let s = solver(&["H2", "O2"]);
        let ox = [0.0, 1.0];
        let fu = [1.0, 0.0];
        let r = s.chem_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let phi = s.weight_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let mut actual = vec![r, phi];
        actual.extend(s.of_ratio_to_weights(&ox, &fu, 2.5).unwrap());
        let amounts = [1.0, 2.0];
        actual.extend(s.weights_to_moles(&amounts).unwrap());
        actual.extend(s.moles_to_weights(&amounts).unwrap());
        for (a, b) in actual.into_iter().zip([10.582243652069236, 10.582243652069236, 9.71848, 24.2962, 0.49606127348850126, 0.06250234383789392, 2.01588, 63.9976]) { close(a, b, 1e-12); }
        close(s.of_ratio_to_chem_eq_ratio(&ox, &fu, r).unwrap(), 0.75, 1e-12);
        close(s.of_ratio_to_weight_eq_ratio(&ox, &fu, phi).unwrap(), 0.75, 1e-12);
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let ox = [0.0, 4.0, 13.0];
        let fu = [1.0, 0.0, 0.0];
        let r = s.chem_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let phi = s.weight_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let mut actual = vec![r, phi];
        actual.extend(s.of_ratio_to_weights(&ox, &fu, 2.5).unwrap());
        let amounts = [1.0, 2.0, 3.0];
        actual.extend(s.weights_to_moles(&amounts).unwrap());
        actual.extend(s.moles_to_weights(&amounts).unwrap());
        for (a, b) in actual.into_iter().zip([22.60582645471247, 22.60582645471247, 145.20338857142858, 85.41375798319328, 277.59471344537815, 0.06233457960936166, 0.06250234383789392, 0.10709160615990919, 16.04246, 63.9976, 84.0402]) { close(a, b, 1e-12); }
        close(s.of_ratio_to_chem_eq_ratio(&ox, &fu, r).unwrap(), 0.75, 1e-12);
        close(s.of_ratio_to_weight_eq_ratio(&ox, &fu, phi).unwrap(), 0.75, 1e-12);
    }
    {
        let s = solver(&["O2(L)", "RP-1"]);
        let ox = [1.0, 0.0];
        let fu = [0.0, 1.0];
        let r = s.chem_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let phi = s.weight_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let mut actual = vec![r, phi];
        actual.extend(s.of_ratio_to_weights(&ox, &fu, 2.5).unwrap());
        let amounts = [1.0, 2.0];
        actual.extend(s.weights_to_moles(&amounts).unwrap());
        actual.extend(s.moles_to_weights(&amounts).unwrap());
        for (a, b) in actual.into_iter().zip([4.5408883574769074, 4.5408883574769074, 32.83927357142857, 13.13570942857143, 0.03125117191894696, 0.14310058762109798, 31.9988, 27.952366]) { close(a, b, 1e-12); }
        close(s.of_ratio_to_chem_eq_ratio(&ox, &fu, r).unwrap(), 0.75, 1e-12);
        close(s.of_ratio_to_weight_eq_ratio(&ox, &fu, phi).unwrap(), 0.75, 1e-12);
    }
    {
        let s = solver(&["N2O4(L)", "N2H4(L)", "C2H8N2(L),UDMH"]);
        let ox = [1.0, 0.0, 0.0];
        let fu = [0.0, 0.5, 0.5];
        let r = s.chem_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let phi = s.weight_eq_ratio_to_of_ratio(&ox, &fu, 0.75).unwrap();
        let mut actual = vec![r, phi];
        actual.extend(s.of_ratio_to_weights(&ox, &fu, 2.5).unwrap());
        let amounts = [1.0, 2.0, 3.0];
        actual.extend(s.weights_to_moles(&amounts).unwrap());
        actual.extend(s.moles_to_weights(&amounts).unwrap());
        for (a, b) in actual.into_iter().zip([2.9984410200619216, 2.9984410200619216, 98.63052857142858, 19.726105714285715, 19.726105714285715, 0.01086826575083414, 0.062411921176240025, 0.0499182007084391, 92.011, 64.09032, 180.29496]) { close(a, b, 1e-12); }
        close(s.of_ratio_to_chem_eq_ratio(&ox, &fu, r).unwrap(), 0.75, 1e-12);
        close(s.of_ratio_to_weight_eq_ratio(&ox, &fu, phi).unwrap(), 0.75, 1e-12);
    }
    // Oxygen within the fuel makes chemical r differ from mass-based phi.
    let s = solver(&["CH3OH", "O2"]);
    close(s.chem_eq_ratio_to_of_ratio(&[0.0,1.0], &[1.0,0.0], 0.75).unwrap(), 2.163754954716528, 1e-12);
    close(s.weight_eq_ratio_to_of_ratio(&[0.0,1.0], &[1.0,0.0], 0.75).unwrap(), 1.9973122658921798, 1e-12);
    let v = s.get_valence().unwrap();
    assert_eq!(v["C"], (0.0,4.0));
    assert_eq!(v["H"], (0.0,1.0));
    assert_eq!(v["O"], (-2.0,0.0));
}
#[test]
fn conversions_validate_streams_and_keep_amount_scale() {
    let s = solver(&["H2", "O2"]);
    for bad in [
        vec![],
        vec![0.0],
        vec![-1.0, 1.0],
        vec![f64::NAN, 1.0],
        vec![f64::INFINITY, 1.0],
    ] {
        assert!(s.weights_to_moles(&bad).is_err());
        assert!(s.moles_to_weights(&bad).is_err());
        assert!(s.of_ratio_to_weights(&bad, &[1.0, 0.0], 1.0).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            s.of_ratio_to_weights(&[0.0, 1.0], &[1.0, 0.0], bad)
                .is_err()
        );
    }
    assert!(
        s.of_ratio_to_weights(&[0.0, 0.0], &[1.0, 0.0], 1.0)
            .is_err()
    );
    assert!(
        s.chem_eq_ratio_to_of_ratio(&[1.0, 0.0], &[1.0, 0.0], 1.0)
            .is_err()
    );
    assert!(
        s.weight_eq_ratio_to_of_ratio(&[0.0, 1.0], &[1.0, 0.0], 0.0)
            .is_err()
    );
    assert_eq!(s.weights_to_moles(&[0.0, 0.0]).unwrap(), vec![0.0, 0.0]);
    let w = [3.0, 7.0];
    for (actual, expected) in s
        .moles_to_weights(&s.weights_to_moles(&w).unwrap())
        .unwrap()
        .into_iter()
        .zip(w)
    {
        close(actual, expected, 1e-15);
    }
    let fuel = s
        .of_ratio_to_weights(&[0.0, 1.0], &[1.0, 0.0], 0.0)
        .unwrap();
    assert!(fuel[0] > 0.0 && fuel[1] == 0.0);
}
#[test]
fn nasa_temperature_ranges() {
    // cea.Reactant(name).get_valid_temperature_range(); bounds in K.
    for (name, expected) in [
        ("O2(L)", (80.17, 100.17)),
        ("RP-1", (288.15, 308.15)),
        ("N2O4(L)", (288.15, 308.15)),
        ("N2H4(L)", (100.0, 800.0)),
        ("C2H8N2(L),UDMH", (288.15, 308.15)),
    ] {
        assert_eq!(
            Reactant::named(name).unwrap().get_valid_temperature_range(),
            expected
        );
    }
    // A reference-only reactant is accepted only inside NASA's nominal range.
    assert!(Reactant::named("RP-1").unwrap().properties(320.0).is_err());
}
#[test]
#[rustfmt::skip]
fn nasa_extra_state_properties() {
    // cea.EqSolver.solve: T in K, P in bar, energy/entropy inputs divided by
    // cea.R. Read gibbs_energy/cv_fr/cv_eq (kJ units), volume, M and MW.
    // HP/SP/TV/UV/SV inputs come from native TP at 2500 K and 1 MPa.
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::TP(1000.0, 100000.0), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-24784151.972973328, 1823.4387563274756, 1823.4681546990873, 4.601065377453583, 18.070836464839783, 18.070836464839783]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::TP(2500.0, 1000000.0), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.42368303, 2561.410028200984, 3967.2804551352174, 1.1658151744039886, 17.829820246272554, 17.829820246272554]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::HP(-7414587.537578467, 1000000.0), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.47674519, 2561.4100290380857, 3967.2804713846836, 1.165815176398601, 17.829820242621423, 17.829820242621423]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::SP(14447.442754441825, 1000000.0), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.681826375, 2561.4100325731474, 3967.2805324947376, 1.1658151838422421, 17.8298202297576, 17.8298202297576]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::TV(2500.0, 1.1658151744039886), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.42377824, 2561.410028180387, 3967.2804550995875, 1.1658151744039886, 17.829820245253917, 17.829820245253917]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::UV(-8580402.711982455, 1.1658151744039886), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.482192695, 2561.41002915989, 3967.280472767975, 1.1658151744039886, 17.829820242143732, 17.829820242143732]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::SV(14447.442754441825, 1.1658151744039886), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-43533194.709011786, 2561.410033216861, 3967.2805392199407, 1.1658151744039886, 17.82982022778365, 17.82982022778365]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["H2", "O2"]);
        let st = s.solve(Mode::TP(3500.0, 5000000.0), &[1.0, 8.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-56211654.147958845, 2704.1720218362093, 8423.950001826235, 0.3589928160364408, 16.21246091846363, 16.21246091846363]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::TP(1000.0, 100000.0), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-10856389.119452825, 1016.7704677771105, 1016.8292317720933, 3.0089421874878965, 27.632667834477783, 27.632667834477783]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::TP(2500.0, 1000000.0), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261977.63757231, 1230.398051706898, 1938.3096222175946, 0.7600504833440307, 27.348545202610264, 27.348545202610264]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::HP(210197.69701224382, 1000000.0), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261977.6902833, 1230.398051968683, 1938.309632391966, 0.7600504852056996, 27.3485451970386, 27.3485451970386]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::SP(9388.870133833823, 1000000.0), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261981.98676664, 1230.3980733058088, 1938.3104617073795, 0.7600506369518756, 27.34854474286008, 27.34854474286008]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::TV(2500.0, 0.7600504833440307), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261977.637572262, 1230.3980517068965, 1938.3096222175975, 0.7600504833440307, 27.348545202610598, 27.348545202610598]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::UV(-549852.7863317868, 0.7600504833440307), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261977.69788883, 1230.3980520224152, 1938.3096336659628, 0.7600504833440307, 27.348545196267853, 27.348545196267853]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::SV(9388.870133833823, 0.7600504833440307), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-23261982.61471058, 1230.3980777532404, 1938.3105668808682, 0.7600504833440307, 27.348544679074724, 27.348544679074724]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CH4", "O2", "N2"]);
        let st = s.solve(Mode::TP(3500.0, 5000000.0), &[1.0, 4.0, 13.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-31350885.509877317, 1250.5654396455843, 3134.1419918551856, 0.22487333852504318, 25.88193441772483, 25.88193441772483]) { close(a,b,1e-4); }
    }
    {
        let s = solver(&["CO"]);
        let st = s.solve(Mode::TP(800.0, 100000.0), &[1.0]).unwrap();
        let actual = [mixture::gibbs_energy(&st), mixture::cv_frozen(&st), mixture::cv_equilibrium(&st).unwrap(), mixture::specific_volume(&st), mixture::molecular_weight_m(&st), mixture::molecular_weight_mw(&st)];
        for (a,b) in actual.into_iter().zip([-10453079.838419281, 1109.4680682064347, 3364.2400590837524, 1.247036344971401, 53.339327492917185, 28.010099999999998]) { close(a,b,1e-4); }
    }
}
