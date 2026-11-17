# Implementation and remaining gaps

The target is **verbatim NASA CEA 3.3.4 behavior**, including retained fields,
failed solves and reported successes that do not satisfy physical checks.
The source reference is commit
[`2f79a647e85737742fb08a127da00d1ea1d069da`](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da).
The native PyPI wheel `cea==3.3.4` is the numerical oracle.
This is still an experimental port; the tests establish parity for their inputs.

| Capability | Implementation and evidence |
|---|---|
| Thermodynamic fits | NASA interval selection, gas extensions, reference enthalpies, species order and polynomial evaluation; bundled NASA data remain unchanged |
| TP, HP, SP, TV, UV, SV | NASA Newton equations and limits; Rust State exposes `converged` and `last_error`, along with the fields NASA retains |
| Condensed products | Compiled phase numbers, activation/removal order, coexistence, singular recovery, pivot and SMNOL limits; no added recovery solves; [details](condensed-phase-insertion.md) |
| Ions | NASA truncation, charge polish, ion shutdown and SMNOL behavior |
| IAC/FAC rockets | Pressure and area schedules, native stops, frozen reference states and partial output; [details](nozzle-continuation.md) |
| Incident/reflected shocks | Analytic NASA volume partials and damped updates; the 1000 m/s RP-1311 case retains its last valid state and code 9 |
| Equilibrium detonation | NASA eight-step iteration; final T and wave speeds update while preceding thermodynamic fields remain stored; nonconverged output remains available |
| Transport | NASA gas-species selection, missing-data estimates, reaction systems and continuation after singular elimination; station frozen Cp is replaced where NASA replaces it |
| Condensed transport | NASA's transport calculation uses the gas block (`num_gas`), so condensed species are excluded from transport Cp and mixing |
| Extra state fields | Stored Cp, Cv, Gibbs energy and M/MW preserve native frozen and retained reporting values |
| Sensitivities | Supported hard-truncation derivatives use NASA's scaled multiple-RHS pivoting; singular systems return zero updates; closure diagnostics do not reject output |
| Input decks | Structured Rust/Python input and case results for all 14 sample decks; parser and report gaps are listed in [input-decks.md](input-decks.md) |
| NASA Python API | Public names and arguments exist under `cea_rs.nasa`; supported methods keep NASA units, shapes and status behavior; [member inventory](nasa-python-api.md) |

## Failure contract

Rust and the original `cea_rs` interface use SI. `cea_rs.nasa` uses NASA's bar
and kJ conventions. Numerical failure returns a state or partial flow result.
Check `State.converged`, `ShockResult.converged`, `DetonationResult.converged`,
`RocketResult.complete`, and the NASA adapter's `converged`/`last_error`.
Codes are 0 (success), 8 (not converged), 9 (last valid shock state) and
10 (Fortran abort). The adapter warns on numerical nonconvergence and retains
the previous report. Fortran aborts raise RuntimeError after updating the
solution, retaining the native converged flag and working iterate. Input validation still
rejects malformed or nonfinite requests.

A native success is not a proof of conservation or phase stability. The port
returns native false positives, including Air/Jet-A at 20 MPa and Br2 at 262 K.
NaNs and stale or zero fields are part of some retained native flow results.
Physical checks can be used by a consumer, but do not change solver acceptance.

## Remaining gaps

These features raise `NotImplementedError` by name in the NASA Python API:

- Smooth truncation and a nondefault truncation width.
- Custom binary databases, `init_thermo` and `init_trans`.
- Overrides of fitted named reactant data.
- Frozen-reactant rockets (`n_frz=-1`) and frozen detonation.
- MATLAB wrappers and native diagnostic logging other than `LOG_NONE`.
- Untruncated frozen rocket logs, and logs lost through element reduction.
- Analytic sensitivities for ions, reduced element systems, absent-element
  weight perturbations, or UV/SV problems with condensed phases.
- Sensitivities of a nonconverged solution and verbose FD diagnostics.
- Report fields on unsolved flow solutions, and uncomputed derivative workspaces.
  The native flow storage is uninitialized and can expose garbage or crash.
- Single-product pools with incompatible element inventories. The reference
  wheel can crash on this input; the safe Rust implementation names the gap.

The port has no NASA CLI, report/plot generator, C/Fortran/Excel ABI or custom
database loader. Deck formatting preferences are metadata. Broader chemistry,
Linux/Windows rounding and production consumers have not been exhaustively
validated. EqSolution guesses are supported, including NASA's setter quirk:
the first solve restores the constructor seed even after T_init/nj_init setters.
NASA's broken shock amount getters still raise AttributeError.

[Current validation](verbatim-parity.md) records commands and measured differences.
[Future improvements](future-improvements.md) lists the removed corrections.
