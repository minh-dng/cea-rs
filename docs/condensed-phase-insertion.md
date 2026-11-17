# Condensed phases: NASA behavior

The target is the logic in NASA `source/equilibrium.f90` at
[2f79a647](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da).
The port keeps its quirks. A physically better equilibrium is a different feature.

`check_condensed_phases`, `test_condensed`, `correct_singular` and `gauss`
define the active-set behavior. The Rust code follows their phase numbering,
activation order, removal order, singular recovery and stopping limits.

- A compiled phase number counts fit intervals. Adjacency uses the difference
  between those numbers. Ice and liquid water therefore need not be adjacent.
- Candidates are scanned in NASA's melting-temperature order. The current phase
  stays in that list. The scan stops when it reaches that phase in range.
- Phase changes transfer the amount or create a coexistence pair and pin its
  temperature, including during fixed-temperature solves. There is no separate
  Gibbs-ranking or fit-range-ranking rule for a sibling replacement.
- Insertion uses the most negative Gibbs decrease per molecular weight. It
  applies NASA's lower/upper bounds, including the special 200 K lower bound.
  A reversed Br2(cr) interval is not repaired.
- Removal uses the lowest nonpositive `nj * Cp`, including NASA's tie order.
  Newly inserted phases go to the front of the active order.
- Gaussian pivot ratios are limited to `1e25`. The final diagonal division is
  performed without an extra zero or finite-value check.
- Singular recovery changes `xsize` and `tsize` to 80, removes the phase NASA
  selects, or reduces element equations in NASA's order. Gas reseeding uses
  `1e-6` and `SMNOL=-13.815511`. There are at most eight singular recoveries.
- The iteration and repeated-convergence limits come from NASA. No cold or
  fixed-temperature recovery solve is added.

Product order and arithmetic also affect phase selection. Rust uses NASA's
species-name order and fused multiply-add operations for the polynomial fits,
Gaussian elimination, matrix products and updates used by the reference build.
For Ar/NaCN at 300 K, rounding a trace graphite amount below zero used to
remove it and change the subsequent solve. The pinned Fortran source was built
in scratch space and stepped through to identify that difference. The port now
follows the native insertion/removal path without a recovery retry.

## Native regressions

`tests/data/nasa/retained-3.3.4.json` records the wheel's values and status codes.
The Python tests use these records with or without an installed native wheel.
When `cea==3.3.4` is installed they also repeat the calls against that wheel.

| Input, with argon | Native behavior retained by the port |
|---|---|
| Hg HP, H=-15551.803948663037 J/kg, 1 bar, weights 1:20 | Error 8; final T and mole amounts remain readable; thermodynamic report fields remain zero |
| Water TP, 273.15 K, 1 bar, weights 1:1 | Error 8; no lower-Gibbs fixed-T recovery |
| NaCN TP, 288.5 K, 1 bar, weights 1:1 | Error 8; transition cycling is retained |
| NaCN TP, 300 K, 1 bar, weights 1:1 | Success with NaCN(III), without a retry |
| Br2 TP, 262 K, 1 bar, weights 1:20 | Success with supersaturated gas and no solid |

The old tests that demanded a corrected cold solve, phase adjacency, Br2 solid
or a fixed-T recovery have been removed. See [future improvements](future-improvements.md)
for the work that can return later. Broader metal systems and cross-platform
rounding have not been exhaustively compared.
