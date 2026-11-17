//! Goldens: cea==3.3.4, NASA source 2f79a647e85737742fb08a127da00d1ea1d069da.
//! Recorded using an independent throwaway uv project (cea==3.3.4, numpy):
//! ```python
//! cea.init(path=str(Path(cea.__file__).parent / "data"), trans=False)
//! r = cea.Mixture(names)
//! p = cea.Mixture(names, products_from_reactants=True, ions=False)
//! s = cea.DetonationSolver(p, reactants=r, transport=False)
//! out = cea.DetonationSolution(s)
//! weights = r.moles_to_weights(np.array(moles))
//! s.solve(out, weights, T1, p1_bar)
//! assert out.converged and not out.last_error
//! ```
//! H2/O2 uses moles [2, 1] at 298.15/500 K and 1/20 bar (RP-1311 example 6).
//! CH4/air uses [phi, 2, 7.52] CH4/O2/N2; acetylene/O2 uses [phi, 2.5].
//! Both use phi=0.5/1/2, 300 K, 1 bar. Full product pools, including condensates.
//! Named Air also uses CH4/Air moles [1, 4/0.41895] at 300 K, 1 bar; the
//! oxygen allowance accounts for the trace carbon in Air's NASA formula.
//! All cases, including named Air, use 1e-4.
//! Native Python rounds T1/p1 to float32; goldens use its recorded effective T1.
//! A separate check covers exactly 298.15 K without that wrapper rounding.
use cea_rs::{Error, detonation::DetonationResult, solver::Solver, thermo::Reactant};

struct Golden {
    names: &'static [&'static str],
    weights: &'static [f64],
    values: [f64; 23],
    tolerance: f64,
}

const FIELDS: [&str; 23] = [
    "P1", "T1", "H1", "M1", "gamma1", "sound1", "T", "P", "density", "H", "U", "S", "Cp_eq",
    "Cp_fr", "gamma_s", "M", "sound", "P/P1", "T/T1", "M/M1", "rho/rho1", "Mach", "velocity",
];

fn values(r: &DetonationResult) -> [f64; 23] {
    let a = &r.unburned;
    let b = &r.burned;
    [
        a.pressure_pa,
        a.temperature_k,
        a.enthalpy_j_kg,
        a.molecular_weight_kg_kmol,
        a.gamma_s,
        a.sound_speed_m_s,
        b.temperature_k,
        b.pressure_pa,
        b.density_kg_m3,
        b.enthalpy_j_kg,
        b.energy_j_kg,
        b.entropy_j_kg_k,
        b.cp_equilibrium_j_kg_k,
        b.cp_frozen_j_kg_k,
        b.gamma_s,
        b.molecular_weight_kg_kmol,
        b.sound_speed(),
        r.pressure_ratio,
        r.temperature_ratio,
        r.molecular_weight_ratio,
        r.density_ratio,
        r.mach,
        r.velocity_m_s,
    ]
}

fn solver(names: &[&str]) -> Solver {
    Solver::new(
        names.iter().map(|s| Reactant::named(s).unwrap()).collect(),
        true,
        false,
        None,
        vec![],
    )
    .unwrap()
}

#[rustfmt::skip]
const GOLDENS: &[Golden] = &[
    Golden {
        tolerance: 1e-4,
        names: &["H2", "O2"],
        weights: &[4.03176, 31.9988],
        values: [
            100000.0, 298.1499938964844, -0.015252560203935608, 12.010186666666668,
            1.4016185359557258, 537.8680777057075, 3674.280834212604, 1876844.7092057355,
            0.8908078132439363, 2831087.377492153, 724185.4382611987, 17425.307874279137,
            16264.788457436047, 3256.7732538500895, 1.1287876732241022, 14.499872669936801,
            1542.1568188645151, 18.76874735586261, 12.323598555860743, 1.207297860755159,
            1.8386781913548433, 5.271794754190256, 2835.5301104953455,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["H2", "O2"],
        weights: &[4.03176, 31.9988],
        values: [
            2000000.0, 298.1499938964844, -0.015252560203935608, 12.010186666666668,
            1.4016185359557258, 537.8680777057075, 4283.381543265734, 41444610.47686056,
            17.752061500409685, 3146394.129969522, 811757.6441786132, 15772.994407035496,
            10269.02784434188, 3320.467079464718, 1.143789245825163, 15.254689005274097,
            1634.1165777994192, 20.722732270250138, 14.366532386221998, 1.2701458710556341,
            1.832063428672189, 5.56605856429226, 2993.805220373268,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["H2", "O2"],
        weights: &[4.03176, 31.9988],
        values: [
            100000.0, 500.0, 495429.1155003234, 12.010186666666668,
            1.3857934682370656, 692.592379019447, 3600.031140025505, 1095661.960986513,
            0.5219118523374723, 3172553.9207331794, 1073230.1439701081, 17830.437511897715,
            18059.168627332307, 3249.0002008466213, 1.1266107198759425, 14.25814551676691,
            1537.895679914176, 10.956873530770935, 7.20006228005101, 1.187171016778555,
            1.8065669734435592, 4.011467102580502, 2778.311543934478,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["H2", "O2"],
        weights: &[4.03176, 31.9988],
        values: [
            2000000.0, 500.0, 495429.1155003234, 12.010186666666668,
            1.3857934682370656, 692.592379019447, 4209.491060086682, 24337988.20630646,
            10.420347682524882, 3500315.678458595, 1164694.1667962046, 16149.066337511922,
            11388.40871192383, 3313.5204746938502, 1.142311542792781, 14.985234989169017,
            1633.4040555167624, 12.169063563721263, 8.418982120173364, 1.2477104149229723,
            1.8034708246936064, 4.253290461023082, 2945.796559060697,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["CH4", "O2", "N2"],
        weights: &[8.02123, 63.9976, 210.660768],
        values: [
            100000.0, 300.0, -130016.64537684346, 28.211536726546907,
            1.392302304394853, 350.8587861788471, 2008.927603374277, 1158075.3066452267,
            1.955307736255358, 608302.4142654123, 16029.75806292841, 8682.43353017193,
            1467.6699046600775, 1381.0435141902508, 1.2536251374776102, 28.201956092684846,
            861.6774010645166, 11.580754854892403, 6.696425344580923, 0.9996604001421502,
            1.7288061140116187, 4.245791235527534, 1489.673159265978,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["CH4", "O2", "N2"],
        weights: &[16.04246, 63.9976, 210.660768],
        values: [
            100000.0, 300.0, -254626.06064705394, 27.633158555133083,
            1.3868783073509923, 353.820401112642, 2776.767255191866, 1706310.0769346985,
            2.0006307920898103, 871806.0382630291, 18919.996679950145, 9465.242449838295,
            2723.9512738818635, 1546.9538801544675, 1.1681231872823734, 27.069804207075503,
            998.136264477763, 17.063105206612114, 9.255890850639553, 0.9796131033325928,
            1.8059026470626216, 5.094496859087938, 1802.5369221495891,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["CH4", "O2", "N2"],
        weights: &[32.08492, 63.9976, 210.660768],
        values: [
            100000.0, 300.0, -484293.9346876957, 26.627021527777778,
            1.3777918129103375, 359.2604938318731, 2160.3578842656475, 1466692.010338513,
            1.8520959233183758, 524802.5336687126, -267106.7793816935, 10400.017851507495,
            1746.3033562225892, 1702.5530375828257, 1.267821775528782, 22.682338116152238,
            1001.9968920722694, 14.666708013124353, 7.201192947552158, 0.8518541246714216,
            1.7349972913033194, 4.839001013157038, 1738.4618936397317,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["C2H2,acetylene", "O2"],
        weights: &[13.01864, 79.997],
        values: [
            100000.0, 300.0, 1228574.6873416915, 31.005213333333334,
            1.3533336998742547, 329.9620883756228, 3789.013285786212, 2598941.470578519,
            2.292718785577601, 2778740.658244432, 1645177.53357026, 9819.928135800421,
            7027.161652048731, 1522.581221382826, 1.1385979061353633, 27.791826445398026,
            1136.0777124112, 25.98940744113611, 12.630044285954039, 0.8963597878399168,
            1.8444801264481476, 6.35064704845026, 2095.4727626431327,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["C2H2,acetylene", "O2"],
        weights: &[26.03728, 79.997],
        values: [
            100000.0, 300.0, 2154184.719369856, 30.29550857142857,
            1.3287958903888124, 330.7645609175263, 4210.5615197099105, 3362867.4258357165,
            2.2368564409206355, 4226752.691315341, 2723362.771284692, 11175.125440204285,
            8867.539497106163, 1690.2211220841355, 1.1527533783237112, 23.286548505565364,
            1316.4488181219792, 33.62860751112237, 14.035205065699703, 0.76864689201906,
            1.8416952997586242, 7.32997995305965, 2424.497600708045,
        ],
    },
    Golden {
        tolerance: 1e-4,
        names: &["C2H2,acetylene", "O2"],
        weights: &[52.07456, 79.997],
        values: [
            100000.0, 300.0, 3457965.9349399633, 29.349235555555556,
            1.300931055491797, 332.5122892999988, 4524.567183692885, 4343918.335946467,
            2.1454762649264216, 6250428.754049346, 4225741.551977018, 12736.763833622466,
            7007.94561074202, 1920.711563556846, 1.1865084533127537, 18.580436787412324,
            1549.9379406614728, 43.43903951635556, 15.081890612309616, 0.6330807748720123,
            1.8234121116094995, 8.499462137759926, 2826.175613245215,
        ],
    },
    Golden {
        // Named-Air equilibrium properties differ even at the same native TP.
        // Its nearly zero U magnifies a 3.78 J/kg difference; all requested cases above use 1e-4.
        tolerance: 1e-4,
        names: &["CH4", "Air"],
        weights: &[16.04246, 276.54962071846285],
        values: [
            100000.0, 300.0, -257075.8793088762, 27.73995003722604,
            1.3881450469533316, 353.2999239888086, 2778.4490521204557, 1707261.9977507244,
            2.008119349888232, 865733.1632854218, 15553.617012065615, 9424.251859887034,
            2710.488869776368, 1535.8471124708508, 1.1684806778974983, 27.17242676692547,
            996.7037745547788, 17.072624464263413, 9.261496840401518, 0.9795413016411718,
            1.8056840470260094, 5.094063098018617, 1799.7321053241722,
        ],
    },
];

#[test]
fn matches_native_cj_states() {
    let mut worst = 0.0_f64;
    for g in GOLDENS {
        let r = solver(g.names)
            .detonation(g.weights, g.values[1], g.values[0], false)
            .unwrap_or_else(|e| {
                panic!("{:?}, {} K, {} Pa: {e}", g.names, g.values[1], g.values[0])
            });
        for (i, actual) in values(&r).into_iter().enumerate() {
            let expected = g.values[i];
            let relative = ((actual - expected) / expected).abs();
            worst = worst.max(relative);
            assert!(
                relative <= g.tolerance,
                "{:?}, {} K, {} Pa, {}: {actual} vs {expected}, relative {relative}",
                g.names,
                g.values[1],
                g.values[0],
                FIELDS[i]
            );
        }
    }
    eprintln!("worst native relative difference: {worst:e}");
}

#[test]
fn exact_reference_temperature_uses_fixed_reactant_fits() {
    // Independent native Mixture.calc_property(cea.ENTHALPY, weights,
    // np.array([298.15, 298.15])) in J/kg, avoiding DetonationSolver's float32 T1.
    let weights = [4.03176, 31.9988];
    let s = solver(&["H2", "O2"]);
    for golden in &GOLDENS[..2] {
        let r = s
            .detonation(&weights, 298.15, golden.values[0], false)
            .unwrap();
        let mut expected = golden.values;
        expected[1] = 298.15;
        expected[2] = -0.0005062607273303999;
        for (i, actual) in values(&r).into_iter().enumerate() {
            assert!(
                ((actual - expected[i]) / expected[i]).abs() < 1e-4,
                "exact 298.15 K, {} Pa, {}: {actual} vs {}",
                golden.values[0],
                FIELDS[i],
                expected[i]
            );
        }
        assert!(r.burned.equilibrium);
    }
}

#[test]
fn rejects_unsupported_frozen_and_invalid_initial_gas() {
    let s = solver(&["H2", "O2"]);
    // Native 3.3.4 solve(..., frozen=True) returns converged=False, last_error=8
    // and warns CEA_NOT_CONVERGED. Its Fortran frozen branch calls abort.
    assert!(matches!(
        s.detonation(&[1.0, 8.0], 298.15, 1e5, true),
        Err(Error::Unsupported(_))
    ));
    for weights in [vec![], vec![0.0, 0.0], vec![-1.0, 8.0], vec![f64::NAN, 1.0]] {
        assert!(matches!(
            s.detonation(&weights, 300.0, 1e5, false),
            Err(Error::Input(_))
        ));
    }
    for (t, p) in [
        (0.0, 1e5),
        (300.0, 0.0),
        (f64::NAN, 1e5),
        (300.0, f64::INFINITY),
    ] {
        assert!(matches!(
            s.detonation(&[1.0, 8.0], t, p, false),
            Err(Error::Domain(_))
        ));
    }
    assert!(matches!(
        solver(&["H2O(L)"]).detonation(&[1.0], 300.0, 1e5, false),
        Err(Error::Input(_))
    ));
    assert!(matches!(
        solver(&["C2H5OH(L)"]).detonation(&[1.0], 300.0, 1e5, false),
        Err(Error::Input(_))
    ));
    let retained = solver(&["N2"])
        .detonation(&[1.0], 300.0, 1e5, false)
        .unwrap();
    assert!(!retained.converged);
    assert_eq!(retained.last_error, 8);
    assert!((retained.mach / 1.0049374234358726 - 1.0).abs() < 1e-4);
}
