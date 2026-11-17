# Upstream source audit

The behavior reference is NASA CEA 3.3.4 at
[`2f79a647e85737742fb08a127da00d1ea1d069da`](https://github.com/nasa/cea/tree/2f79a647e85737742fb08a127da00d1ea1d069da),
under Apache-2.0. Source was fetched and built in scratch space because the
external drive configured in mise was unavailable.

The audit covered equilibrium phase insertion/removal, compiled phase numbers,
coexistence, Gaussian pivots and singular recovery; rocket, shock and detonation
iterations; transport selection, estimates and Cp replacement; and Python names,
signatures, units, status codes and retained fields. Source and documentation
searches included `Unlike NASA`, `instead of`, `deviation` and retry paths.
Corrections that changed NASA's answer or status were removed and listed in
[future improvements](future-improvements.md).

The prior Rust sources remain attributed in [provenance](provenance.md):
Bit Mango's cea-rs `914ffd48f3acd49dbacfdc54e9e12bf848a5e9ad` and
thermo-species `d4061dddda7f04e2edf75cdb235d19d1fed4161b`, used under
Apache-2.0 with original licenses retained. aselimov/cea-rs was inspected but
not copied. All 11 NASA condensed records with empty or reversed intervals
are included; their intervals are not repaired.

The native wheel supplies independent comparisons and snapshots. Its crashes
on incompatible single-product inventories and uninitialized flow storage are
named unsupported features. [Current validation](verbatim-parity.md) records
what ran; the [implementation matrix](implementation-matrix.md) lists remaining gaps.
