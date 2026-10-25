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

A final focused free-model review is required after these corrections and
before release; its result is recorded below when completed.
