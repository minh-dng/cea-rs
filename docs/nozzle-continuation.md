# Nozzle continuation and NASA stops

The solver follows the pinned NASA 3.3.4 rocket and equilibrium routines.
It carries temperature, gas logs, amounts and active phases between stations.
A failed equilibrium solve restores the saved stable starting point for the
next solve while keeping the report values NASA retains.

Throat calculations use NASA's 22-iteration limit, `4e-5` velocity stopping
criterion and pressure update. Assigned areas use ten iterations and stop on
relative area error **or** a small pressure update, both at `4e-5`. There is
no extra root bracket or independent mass-flow acceptance check. Pressure
schedules use NASA's sonic-speed area basis; area schedules retain the throat's
velocity basis. These bases can differ at the native stopping point.

Frozen expansion uses NASA's eight iterations, temperature updates, fit limits
and starting-station selection. Frozen Cv retains the original station's field
or its zero workspace value. FAC chamber closure and the later frozen throat
are separate steps, as in NASA; the port does not reconcile them again.

## Values that must remain inconsistent

- RP-1311 example 12 requests area 75. NASA stops at about 74.99733854. Rust
  keeps that stop and the corresponding composition.
- Air/Jet-A(L), weights 25:1, 20 MPa, area 80 reports success. Its exit retains
  earlier thermodynamic fields while the flow calculation uses later iterate
  amounts. Native energy and mass-flow inconsistencies are reproduced.
- RP-1311 example 9 with `n_frz=3` stops after seven stations with error 8.
  NASA's zero gamma workspace initializes the frozen throat at infinity
  pressure. Its Mach, impulse and characteristic velocity include NaNs.
- Rocket performance recomputes the chamber enthalpy from its stored
  composition. Isp can therefore differ from the velocity used during the
  flow iteration, especially when the energy difference is close to zero.
- FAC contraction ratio one still runs the combustor area iteration. Its
  returned Mach need not equal one.
- Large subsonic-area requests can meet NASA's pressure-update stop before
  reaching their assigned area. The port returns the native result.

The native wheel and stored native snapshots supply the regressions. Corrected
phase, geometry or conservation results belong in
[future opt-in work](future-improvements.md). Current validation and gaps are
listed in [the implementation matrix](implementation-matrix.md).
