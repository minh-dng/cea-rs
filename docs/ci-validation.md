# CI validation layer

The GitLab pipeline runs Rust tests and strict Clippy, then separately builds
and tests a release Python wheel from outside the source directory. Official
Rust 1.99.1 and Python 3.15.1 images are pinned by manifest digest; maturin is
pinned to 1.15.0. This pipeline does not publish packages or deploy consumers.

GitLab CI lint accepted this configuration with no errors or warnings and
returned exactly `rust-checks` and `python-wheel` (2026-10-25).

Both pinned images started locally on Linux aarch64 and reported their expected
versions. The Docker Desktop bind mount of the T7 drive stalled before a
container process began; a separate 15-second mount probe also timed out.
Therefore local Linux project tests and wheel compilation were not completed.
No project build artifacts were moved to the internal disk to work around it.
The existing macOS Rust, wheel and scientific checks remain separate evidence.

Actual GitLab pipeline execution and the independent free CI review are recorded
in the merge request. Lint success is not a successful pipeline run.
