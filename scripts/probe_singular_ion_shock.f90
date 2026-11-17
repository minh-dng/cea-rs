! Native-only rounding probe for the singular-ion-shock-roundoff parity gap.
program probe_singular_ion_shock
    use cea
    use cea_equilibrium, only: dp
    use fb_utils
    implicit none
    type(ThermoDB) :: db
    type(TransportDB) :: tr
    type(Mixture) :: reactants, products
    type(ShockSolver) :: solver
    type(ShockSolution) :: solution
    character(1024) :: data_dir
    real(dp) :: weights(2)

    call get_command_argument(1, data_dir)
    if (len_trim(data_dir) == 0) error stop 'Pass the NASA build directory containing thermo.lib and trans.lib'
    call set_log_level(log_levels%none)
    db = read_thermo(trim(data_dir)//'/thermo.lib')
    tr = read_transport(trim(data_dir)//'/trans.lib')
    reactants = Mixture(db, [character(15) :: 'N2', 'CH4'], ions=.true.)
    products = Mixture(db, reactant_names=[character(15) :: 'N2', 'CH4'], ions=.true.)
    solver = ShockSolver(products, reactants=reactants, trace=1.d-10, ions=.true., all_transport=tr)
    weights = [.985d0, .015d0]*[28.0134d0, 16.04246d0]
    solution = solver%solve(weights, 291.d0, .00029864211683161557d0, &
                            u1=1200.d0, reflected=.false., incident_frozen=.false.)
    if (solution%converged .or. solution%solve_status /= 2) error stop 'Expected a last-valid native shock'
    if (.not. (solution%eq_soln(2)%T > 0.d0 .and. solution%eq_soln(2)%n > 0.d0)) &
        error stop 'Expected positive retained native state'
    print *, 'P21, T21, T2, n2:', solution%p21, solution%t21, solution%eq_soln(2)%T, solution%eq_soln(2)%n
    print *, 'rho2, M2, gamma2:', solution%eq_soln(2)%density, solution%eq_soln(2)%M, solution%eq_soln(2)%gamma_s
    print *, 'u2, sonic2, Mach2:', solution%u(2), solution%v_sonic(2), solution%mach(2)
end program
