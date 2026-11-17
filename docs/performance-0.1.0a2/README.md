# Historical alpha2 performance evidence

These artifacts describe released version 0.1.0a2 and reviewed source snapshot
`4734f3d673d9ccabaf9329099762b4818b849fb7`. They are not measurements of the
current NASA parity implementation. Consumer dependencies were unchanged.

The paired benchmark alternated alpha1 and alpha2 in separate processes for
three rounds, using five warmups and 50 timed warm solves per case, plus fresh
process samples. The raw results remain here. Across nine rocket cases,
median paired speedup was 3.16–12.57 times over alpha1. Alpha2 remained a
median of 7.06–22.37 times slower than native CEA across those rounds.
OS cache and CPU thermal/scheduling state were not controlled.

The bracketed root search and additional physical acceptance checks measured
here have been removed to match NASA CEA 3.3.4. They are listed in
[future improvements](../future-improvements.md). See [current validation](../verbatim-parity.md)
for parity evidence and remaining gaps. No new performance benchmark ran.
