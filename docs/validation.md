# Validation

The current target is verbatim NASA CEA 3.3.4 behavior, including numerical
failures, retained fields and false positives. See [current validation](verbatim-parity.md)
for commands, results, measured differences and remaining gaps.

The 146-row comparison runner isolates native calls in worker processes and
records the loaded native extension and database hashes. Expected values come
from the independent `cea==3.3.4` installation. Rust never reads those values
while solving. [Current comparison data](native-comparison.json).

Older performance and CI artifacts describe their recorded source snapshots.
They do not validate the current numerical implementation. The bracketed roots,
extra recovery solves and stricter acceptance checks from those versions have
been removed; [future improvements](future-improvements.md) records them.
No production consumer or full trajectory acceptance test ran in this pass.
