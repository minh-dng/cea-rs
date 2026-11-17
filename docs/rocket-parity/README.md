# Rocket parity

The current target is the pinned NASA CEA 3.3.4 implementation, including its
stopping rules and retained invalid fields. [Nozzle continuation](../nozzle-continuation.md)
explains the flow rules. [The implementation matrix](../implementation-matrix.md)
lists the remaining unsupported features.

`tests/rocket-native.txt` contains native station goldens. `tests/rocket.rs`,
`tests/nasa_decks.rs` and `tests/test_nasa_api.py` check IAC/FAC schedules,
freezing, transport and partial stops. `retained-3.3.4.json` adds false positives,
failures and NaN workspaces. Stored examples are supplemented by live wheel
comparisons when `cea==3.3.4` is installed.

Throats are not refined past NASA's stopping test. Large subsonic areas and
example 12 retain NASA's small-update stop. Frozen FAC results are not repaired
by an additional combustor/throat closure loop. Transport failures use NASA's
estimates and overwritten matrix values; station frozen Cp follows transport.

`comparison.json` records the current 24-case native sweep. The benchmark,
timing and `existing-comparison-summary.json` files describe earlier source
snapshots; they are historical measurements. Current results and commands
are in [verbatim-parity.md](../verbatim-parity.md).
