# Alpha2 nozzle-root performance evidence

The released 0.1.0a1 artifacts remain unchanged. The reviewed nozzle-root
change is released as version 0.1.0a2. No consumer dependency or production default changes.

A sampled alpha1 NTO/Aerozine50 solve spends substantial time in repeated
thermochemistry calls from bracketed roots. The private root search now uses
secant interpolation with Illinois stationary-endpoint damping and a midpoint
fallback. It preserves the bracket, 80-iteration cap, root tolerances, domain
checks and final sonic/mass/energy closure. The thermochemical core and data
are unchanged. Nonfinite bracket values now fail explicitly.

Three root tests check endpoint/failure behavior, a smooth pressure-like root
with bounded evaluation count, and a strongly skewed polynomial. A nozzle
integration test adds area ratios 1, 1.001 and 100 for all three flow modes.
All existing checks and 146 isolated native comparison rows pass at the
existing limits; maximum rocket station difference remains about 5.84e-5.
Installed-wheel interleaving/thread/closure checks pass on Python 3.12.14 and
3.15.1 outside source.

The paired benchmark alternates released alpha1 and candidate in separate
processes for three rounds, using the same interpreter and inputs, five
warmups and 50 timed warm solves per case, plus fresh-process samples. It
retains all raw results. Across the nine rocket cases, median paired speedup
is 3.16–12.57 times relative to alpha1. The unchanged TP path gives approximately
1.01–1.05 times, consistent with timing noise. Candidate rocket solves remain a median of
7.06–22.37 times slower than native CEA across the three rounds. This is an
improvement to the Rust port; a speed advantage over native CEA is unproven.
OS cache and CPU thermal/scheduling state are not controlled.

Independent free Poolside review completed with scoped approval and no
blockers, against exported commit `4734f3d673d9ccabaf9329099762b4818b849fb7`
without Git metadata. The only clarity note was to label the native timing
ratio as a per-case median. No numerical source changed after that snapshot.
An earlier review was cancelled because it checked out alpha1's rocket file
into the shared worktree; that patch was preserved, the committed candidate
restored, and its rebuilt extension proved byte-identical to the tested wheel.
The cancelled review is not counted as approval.
