# Nozzle continuation checks

The nozzle keeps the previous verified temperature, gas log amounts, stored
amounts, total gas moles, condensed amounts, active phase order and coexistence
pair. NASA resets Newton history and the main trace cutoff on each solve;
the relaxed trace and ion thresholds survive singular recovery. The private
core entry point follows these rules and checks carried phases before assembly.
Public Rust and Python solve signatures are unchanged.

Reference: NASA CEA commit `2f79a647e85737742fb08a127da00d1ea1d069da`,
`source/rocket.f90:RocketSolver_set_init_state` and
`source/equilibrium.f90:EqSolver_solve`. The source was read in a separate clone;
no scripts from that clone were executed. Numerical checks use PyPI `cea==3.3.4`
and its packaged data, independently of Rust. Test comments record the native
API calls and inputs for each golden.

- Air/Jet-A(L), O/F=25, reactants at 298.15 K, 0.7 MPa, area=200:
  chamber and exit match the native rocket within 1e-5 relative. Native stops
  its throat at Mach=0.9999816; Rust preserves its existing sonic closure.
  The raw native throat pressure differs by 2.08e-5 relative and Isp by
  1.62e-5. Refining the native throat with EqSolver SP and Eq. 6.17 until
  `abs(Mach²-1)<1e-10` gives agreement within 1e-5 at every station.
- The same mixture at 20 MPa, area=80 still fails. Warm starts alone do not
  fix the water-phase cycle. Native RocketSolver logs three condensed-set
  convergence warnings but reports success, T=273.412284434 K and
  P=9996.45033548 Pa at the exit. A standalone native SP solve using its
  chamber entropy and exit pressure returns `CEA_NOT_CONVERGED` (error 8).
  Recomputing entropy with the reported composition gives 7457.46148 J/(kg K),
  versus the assigned 7418.56142. Matching that unverified result would violate
  the requirement to return only verified states, so Rust returns an error.
  Ice has compiled phase number 1 and liquid water 3; NASA's phase-number
  test keeps replacing them rather than establishing coexistence.
- RP-1311 example 13 uses N2H4(L)/Be(a)/H2O2(L), weights
  [0.536,0.134,0.33], 298.15 K, 3000 psia, insert BeO(L), trace=1e-10,
  and pressure ratios [3,10,30,300]. All six stations match native at 1e-5
  for T, P, Isp, area, gamma, equilibrium Cp, density and vacuum Isp.
  The explicit trace setting also lets the original checkout solve this deck;
  the older failure report omitted it. Alumina coexistence regressions remain.

The native throat stopping rule and the unverified 20 MPa exit prevent the
requested direct 1e-5 comparison against both unmodified native rocket outputs.
No existing test tolerance was loosened to hide these differences.

A release-build timing on macOS arm64 averaged 100 H2/O2 equilibrium rockets
at weights [1,6], 298.15 K, 3 MPa and area=80. Mean runtime fell from
494.9 µs before continuation to 380.3 µs after it, about 23% faster. This is a
local timing, not a general performance claim.
