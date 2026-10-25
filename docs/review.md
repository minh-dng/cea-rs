# Independent free-model review

The independent implementation review completed through OpenCode using
`opencode/glm-5.3-flash-free`, default options, on 2026-10-25. It inspected Rust
thermochemistry, equilibrium, nozzle closure, Python packaging and the 146-row
native comparison. Task: `cea-rust-implementation-free-review-r1-20261025`.
Acceptance of a task was not treated as completion.

Resolved findings:

- B1: the report could omit the actual native binary. It now identifies and
  hashes the loaded `cea.lib.libcea` extension and packaged NASA binary data,
  rejects a missing extension/data set, and verifies the installed wheel's
  NASA text-data SHA against the source bytes.
- B2: the CEA coefficient format reserves a column between the final polynomial
  and integration constants. The bundled-data parser now rejects nonzero
  content there. This parser supports CEA's NASA-nine-coefficient format; it
  does not accept arbitrary NASA-seven files.
- H2/H3: documented condensed-mole reporting convention and Rust frozen-Cp
  validity; Python exposes equilibrium Cp as None for frozen composition.
- M1: constant-volume entropy convergence now uses pressure recomputed from
  post-update gas moles and temperature. All six-mode comparisons are rerun.
- P1: Cargo.lock committed; exact mise versions and external build locations
  documented. Stale upstream comments claiming TV/UV/SV missing corrected.

Finding H1 is inapplicable to the exposed API: `equilibrium_core` is a private
module (`src/lib.rs`), so its upstream compatibility helper is not accessible
from the rlib or Python. The public Solver checks charge/ion settings before
calling the core. Invalid input/pool tests verify explicit errors. No native
fallback or consumer compatibility claim was added.

The review also noted bracketed rocket-solver cost and exFAT sidecars. These
are retained as documented limitations. The benchmark does not claim speedup.
The unrelated MiMo research task is not counted as this independent review.

Final R3 review completed through OpenCode using
`openrouter/poolside/laguna-s-2.1:free`, default options, with explicit
**scoped alpha approval** and no release blockers. It verified the corrections,
146 passing comparison rows, exact NASA data identity, licenses and source
invariants. No numerical source changed after the first review's corrected
base `73c3808`; later commits preserve data bytes and exclude archive sidecars.
R2 on GLM failed before review with a free-usage-limit error; it was not counted
as a completed review. The R3 task was
`cea-rust-implementation-free-review-r3-20261025`.

Three non-blocking observations were recorded. The comparison oracle runs
Python 3.15.0, while fresh-wheel tests use separate Python 3.12.14 and 3.15.1
interpreters; validation docs now make that distinction explicit. The
`CEA_CARGO_HOME` path label in mise does not configure Cargo by itself; the
published build commands set standard `CARGO_HOME` after mise tool resolution.
Ruff's EXE002 exception is passed explicitly on the command line for exFAT.
These observations do not change the numerical implementation or release scope.
