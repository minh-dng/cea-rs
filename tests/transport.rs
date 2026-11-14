//! Native NASA CEA 3.3.4 goldens, generated 2026-11-15 with trans=True.
//!
//! Generation: uv init --bare /tmp/transport-oracle; uv add cea==3.3.4 numpy.
//! cea.init(path=str(Path(cea.__file__).parent / "data"), trans=True).
//! For each case below construct cea.Mixture(reactants) and either
//! cea.Mixture(only) or cea.Mixture(reactants, products_from_reactants=True),
//! passing ions=True to both mixtures and the solver for ion cases.
//! EqSolver(products, reactants=reactants, transport=True).solve(solution,
//! cea.TP/TV, temperature, pressure_pa/1e5 or volume_m3_kg, numpy.array(weights)).
//! Record viscosity*1e-4, conductivity_fr*0.1, conductivity_eq*0.1, Pr_fr, Pr_eq.
//! The relative tolerance is 1e-4 on every output, with no absolute floor.
//!
//! RP2 uses the sample's explicit product list and stoichiometric H2/Air ratio
//! from Mixture.weight_eq_ratio_to_of_ratio([0,1], [1,0], 1).
//! RP6 goldens come directly from DetonationSolver(..., transport=True) at
//! T1=298.15/500 K, P1=1/20 bar and stoichiometric H2/O2. The Rust solves use
//! the resulting burned T/P printed below, not a Rust detonation calculation.
//! RP11 comes from RocketSolver(..., ions=True, transport=True), Li(cr)/F2(L)
//! mole weights [1,0.5556], inlet T=[298.15,85.02], Pc=1000 psia, iac=True,
//! hc=reactants.calc_property(cea.ENTHALPY, weights, inlet_T)/cea.R, supar=10.
//! Compare its chamber transport at the native T/P using Li/F2 gas reactants
//! with identical element weights. The requested Rust HP solve is blocked by
//! missing named Li(cr) and nonconvergence even with an equivalent custom Li.
use cea_rs::{
    Error,
    solver::{Mode, Solver},
    thermo::Reactant,
};

struct Golden {
    label: &'static str,
    reactants: &'static [&'static str],
    weights: &'static [f64],
    mode: Mode,
    only: Option<&'static [&'static str]>,
    ions: bool,
    expected: [f64; 5],
}

#[test]
fn nasa_transport_goldens() {
    let cases = [
        Golden {
            label: "edge CO 800.0",
            reactants: &["CO"],
            weights: &[1.0],
            mode: Mode::TP(800.0, 100000.0),
            only: None,
            ions: false,
            expected: [
                3.544191322194938e-05,
                0.054973046854013996,
                0.054973046854013996,
                0.7522429798009689,
                0.7522429798009689,
            ],
        },
        Golden {
            label: "edge H2 10000.0",
            reactants: &["H2"],
            weights: &[1.0],
            mode: Mode::TP(10000.0, 100000.0),
            only: Some(&["H", "H2", "H+", "e-"]),
            ions: true,
            expected: [
                9.74715037187639e-05,
                3.1409988728098854,
                6.001023825937937,
                0.6547343141483215,
                0.8264937199244105,
            ],
        },
        Golden {
            label: "edge Ar 300.0",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(300.0, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                2.2730026552534796e-05,
                0.01779377756234932,
                0.01779377756234932,
                0.6646812276034605,
                0.6646812276034605,
            ],
        },
        Golden {
            label: "edge Ar 1000.0",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(1000.0, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                5.5515512393241666e-05,
                0.0434845806180298,
                0.0434845806180298,
                0.664294557174869,
                0.664294557174869,
            ],
        },
        Golden {
            label: "edge Ar 1000.001",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(1000.001, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                5.5515547860089014e-05,
                0.043484609603907404,
                0.043484609603907404,
                0.6642945387442682,
                0.6642945387442682,
            ],
        },
        Golden {
            label: "edge Ar 5000.0",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(5000.0, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                0.00016458268437834122,
                0.12884794055057827,
                0.12884794055057827,
                0.6646427648394841,
                0.6646427648394841,
            ],
        },
        Golden {
            label: "edge Ar 5000.001",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(5000.001, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                0.000164586386080438,
                0.12884868621149315,
                0.12884868621149315,
                0.664653867170347,
                0.664653867170347,
            ],
        },
        Golden {
            label: "edge Ar 16000.0",
            reactants: &["Ar"],
            weights: &[1.0],
            mode: Mode::TP(16000.0, 100000.0),
            only: Some(&["Ar"]),
            ions: false,
            expected: [
                0.00037685935448932435,
                0.2947121942401786,
                0.2947121942401786,
                0.8514478740346828,
                0.8514478740346828,
            ],
        },
        Golden {
            label: "TP H2 500.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                1.7355700813819156e-05,
                0.03595443243922582,
                0.03595443243922582,
                0.9404980300215177,
                0.9404980300215177,
            ],
        },
        Golden {
            label: "TP H2 1000.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(1000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                3.7702140032494544e-05,
                0.0971591336492485,
                0.0971591336492485,
                0.886118911644843,
                0.886118911644843,
            ],
        },
        Golden {
            label: "TP H2 1500.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(1500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                5.552150512039392e-05,
                0.17178089913217737,
                0.1725291678686658,
                0.8455642135989664,
                0.8439545901438994,
            ],
        },
        Golden {
            label: "TP H2 2000.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(2000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                7.117035474267462e-05,
                0.2487192334184982,
                0.28826941020892866,
                0.8186770122867849,
                0.7589796116537159,
            ],
        },
        Golden {
            label: "TP H2 2500.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(2500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                8.523777674276743e-05,
                0.3246460984551095,
                0.6016582853866351,
                0.7949502375999362,
                0.6457109088909009,
            ],
        },
        Golden {
            label: "TP H2 3000.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(3000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                9.805594387005774e-05,
                0.4048219896979334,
                1.5275955759167492,
                0.7604625803726279,
                0.5407024328403954,
            ],
        },
        Golden {
            label: "TP H2 3500.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(3500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00011014301676354204,
                0.5099987330064212,
                3.7395392545364063,
                0.6955314795631424,
                0.4836878187523422,
            ],
        },
        Golden {
            label: "TP H2 4000.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(4000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00012174055919947837,
                0.6677907446633747,
                6.898013939713589,
                0.6014045979428229,
                0.4908428429881078,
            ],
        },
        Golden {
            label: "TP H2 4500.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(4500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00012924709964589327,
                0.8202628748934624,
                7.885335632084871,
                0.5337718161995049,
                0.5211654849000871,
            ],
        },
        Golden {
            label: "TP H2 5000.0",
            reactants: &["H2", "O2"],
            weights: &[1.0, 8.0],
            mode: Mode::TP(5000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00013188398991124874,
                0.897848044911274,
                5.550964542453689,
                0.5082306726107609,
                0.5240184969046531,
            ],
        },
        Golden {
            label: "TP CH4 500.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                2.648876006253898e-05,
                0.04048843156238339,
                0.04167178119184341,
                0.7552679331375964,
                0.7409115757022682,
            ],
        },
        Golden {
            label: "TP CH4 1000.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(1000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                4.3890947312044226e-05,
                0.07584711200592131,
                0.07596484011938447,
                0.7573902637749255,
                0.7572201212343812,
            ],
        },
        Golden {
            label: "TP CH4 1500.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(1500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                5.790549761652158e-05,
                0.11057119286400081,
                0.11068318310486487,
                0.7428550639229904,
                0.7426064377817755,
            ],
        },
        Golden {
            label: "TP CH4 2000.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(2000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                7.042005244071062e-05,
                0.14415165967336893,
                0.15830475214835624,
                0.7250225196381892,
                0.7154467522372568,
            ],
        },
        Golden {
            label: "TP CH4 2500.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(2500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                8.192316731839736e-05,
                0.1760327885086002,
                0.29648670363897683,
                0.7087237169053344,
                0.6297401679027864,
            ],
        },
        Golden {
            label: "TP CH4 3000.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(3000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                9.25317128797211e-05,
                0.20805841474564055,
                0.6690745358743136,
                0.686709924610638,
                0.4835748473347897,
            ],
        },
        Golden {
            label: "TP CH4 3500.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(3500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00010309768326803431,
                0.25010457574632994,
                1.3762698600813097,
                0.6425743518842814,
                0.3740337166926752,
            ],
        },
        Golden {
            label: "TP CH4 4000.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(4000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00011468227207545128,
                0.31207677556949953,
                1.8984798073793763,
                0.5790358792979038,
                0.3603871117736908,
            ],
        },
        Golden {
            label: "TP CH4 4500.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(4500.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00012697331060903475,
                0.37998765889064634,
                1.5822829597816732,
                0.5324439572541297,
                0.4040394521178006,
            ],
        },
        Golden {
            label: "TP CH4 5000.0",
            reactants: &["CH4", "Air"],
            weights: &[1.0, 17.16],
            mode: Mode::TP(5000.0, 1000000.0),
            only: None,
            ions: false,
            expected: [
                0.00013905991441697256,
                0.4362143497149408,
                1.1034301554537145,
                0.5127731135306954,
                0.4625135076045509,
            ],
        },
        Golden {
            label: "RP2 rho=0.091864",
            reactants: &["H2", "Air"],
            weights: &[1.0, 34.29622553441931],
            mode: Mode::TV(3000.0, 10.885657058260037),
            only: Some(&[
                "Ar", "C", "CO", "CO2", "H", "H2", "H2O", "HNO", "HO2", "HNO2", "HNO3", "N", "NH",
                "NO", "N2", "N2O3", "O", "O2", "OH", "O3",
            ]),
            ions: false,
            expected: [
                9.357604933624696e-05,
                0.26316500348995897,
                1.8497787020755163,
                0.6322397677865036,
                0.3558398387784543,
            ],
        },
        Golden {
            label: "RP2 rho=0.0080877",
            reactants: &["H2", "Air"],
            weights: &[1.0, 34.29622553441931],
            mode: Mode::TV(3000.0, 123.64454665726969),
            only: Some(&[
                "Ar", "C", "CO", "CO2", "H", "H2", "H2O", "HNO", "HO2", "HNO2", "HNO3", "N", "NH",
                "NO", "N2", "N2O3", "O", "O2", "OH", "O3",
            ]),
            ions: false,
            expected: [
                9.401386247626087e-05,
                0.3039797432662127,
                4.0256941702095945,
                0.5541276653504491,
                0.3360802412413596,
            ],
        },
        Golden {
            label: "RP2 rho=0.00066054",
            reactants: &["H2", "Air"],
            weights: &[1.0, 34.29622553441931],
            mode: Mode::TV(3000.0, 1513.9128591758258),
            only: Some(&[
                "Ar", "C", "CO", "CO2", "H", "H2", "H2O", "HNO", "HO2", "HNO2", "HNO3", "N", "NH",
                "NO", "N2", "N2O3", "O", "O2", "OH", "O3",
            ]),
            ions: false,
            expected: [
                9.481954200022498e-05,
                0.36153956968438555,
                3.6938323167782254,
                0.47934582639272527,
                0.3990579967002856,
            ],
        },
        Golden {
            label: "RP6 T1=298.15 P1=1.0",
            reactants: &["H2", "O2"],
            weights: &[3.8061863823107354, 30.20849361768926],
            mode: Mode::TP(3674.2808344247833, 1876844.7094525413),
            only: None,
            ions: false,
            expected: [
                0.00011411609913041688,
                0.5398793761233498,
                3.808722117359638,
                0.6883951864638193,
                0.48732187759108153,
            ],
        },
        Golden {
            label: "RP6 T1=298.15 P1=20.0",
            reactants: &["H2", "O2"],
            weights: &[3.8061863823107354, 30.20849361768926],
            mode: Mode::TP(4283.381543661223, 41444610.48358696),
            only: None,
            ions: false,
            expected: [
                0.0001274014898186918,
                0.5945719826677058,
                2.491802212537287,
                0.7114907299065167,
                0.5250374381138327,
            ],
        },
        Golden {
            label: "RP6 T1=500.0 P1=1.0",
            reactants: &["H2", "O2"],
            weights: &[3.8061863823107354, 30.20849361768926],
            mode: Mode::TP(3600.031140211338, 1095661.9611271028),
            only: None,
            ions: false,
            expected: [
                0.00011244301667558024,
                0.5365148935322691,
                4.21081130415499,
                0.68092696622829,
                0.48224129603965954,
            ],
        },
        Golden {
            label: "RP6 T1=500.0 P1=20.0",
            reactants: &["H2", "O2"],
            weights: &[3.8061863823107354, 30.20849361768926],
            mode: Mode::TP(4209.4910604393635, 24337988.210143145),
            only: None,
            ions: false,
            expected: [
                0.00012588666037290883,
                0.5928834027596966,
                2.78640659643428,
                0.7035582810903827,
                0.5145152704848427,
            ],
        },
        Golden {
            label: "RP11 TP chamber",
            reactants: &["Li", "F2"],
            weights: &[6.941, 21.111025635839997],
            mode: Mode::TP(5680.807364834442, 6894757.293168573),
            only: None,
            ions: true,
            expected: [
                0.00014408612213323276,
                0.3151853919141201,
                1.4439192358864181,
                0.7704683680308858,
                0.6727542670630029,
            ],
        },
    ];
    for case in cases {
        let solver = Solver::new(
            case.reactants
                .iter()
                .map(|name| Reactant::named(name).unwrap())
                .collect(),
            true,
            case.ions,
            case.only
                .map(|names| names.iter().map(|s| (*s).into()).collect()),
            vec![],
        )
        .unwrap();
        let state = solver.solve(case.mode, case.weights).unwrap();
        let result = state.transport().unwrap();
        let actual = [
            result.viscosity_pa_s,
            result.conductivity_frozen_w_m_k,
            result.conductivity_equilibrium_w_m_k,
            result.prandtl_frozen,
            result.prandtl_equilibrium,
        ];
        for (i, (actual, expected)) in actual.into_iter().zip(case.expected).enumerate() {
            let relative = (actual / expected - 1.0).abs();
            assert!(
                relative < 1e-4,
                "{} field {i}: {actual} vs {expected}, relative={relative}",
                case.label
            );
        }
    }
}

#[test]
fn transport_rejects_invalid_temperature() {
    let solver = Solver::new(
        vec![Reactant::named("Ar").unwrap()],
        false,
        false,
        None,
        vec![],
    )
    .unwrap();
    let mut state = solver.solve(Mode::TP(1000.0, 1e5), &[1.0]).unwrap();
    state.temperature_k = f64::NAN;
    assert!(matches!(state.transport(), Err(Error::Input(_))));
}

#[test]
fn frozen_nozzle_transport_matches_native() {
    // Same native initialization and SI conversion as above. RocketSolver:
    // H2/O2 weights [1,6], inlet [298.15,298.15] K, Pc=30 bar, supar=10,
    // hc=reactants.calc_property(cea.ENTHALPY, weights, inlet_T)/cea.R.
    // These are the throat and exit columns for n_frz=1 and n_frz=2.
    use cea_rs::rocket::FlowMode;
    let solver = Solver::new(
        vec![
            Reactant::named("H2").unwrap(),
            Reactant::named("O2").unwrap(),
        ],
        true,
        false,
        None,
        vec![],
    )
    .unwrap();
    for (flow, expected) in [
        (
            FlowMode::FrozenAtChamber,
            [
                [
                    9.871962119771526e-5,
                    0.5270891610167598,
                    2.550396496842737,
                    0.696987064065459,
                    0.5074713810713093,
                ],
                [
                    5.662258323536411e-5,
                    0.2612928043839015,
                    4.870970291738992,
                    0.6841900852517316,
                    0.48222226948129954,
                ],
            ],
        ),
        (
            FlowMode::FrozenAtThroat,
            [
                [
                    0.00010202513647301411,
                    0.5447578410440695,
                    2.141170834508275,
                    0.7021095786846847,
                    0.5135800377000277,
                ],
                [
                    5.948960594598093e-5,
                    0.27530253793876164,
                    3.8662328241791166,
                    0.6931034008416365,
                    0.4789807265563807,
                ],
            ],
        ),
    ] {
        let rocket = solver
            .rocket(&[1.0, 6.0], &[298.15, 298.15], 3e6, 10.0, flow)
            .unwrap();
        for (station, expected) in [rocket.throat, rocket.exit].iter().zip(expected) {
            let result = station.state.transport().unwrap();
            for (actual, expected) in [
                result.viscosity_pa_s,
                result.conductivity_frozen_w_m_k,
                result.conductivity_equilibrium_w_m_k,
                result.prandtl_frozen,
                result.prandtl_equilibrium,
            ]
            .into_iter()
            .zip(expected)
            {
                assert!(
                    (actual / expected - 1.0).abs() < 1e-4,
                    "{actual} vs {expected}"
                );
            }
        }
    }
}
