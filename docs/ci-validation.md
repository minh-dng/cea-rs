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

## Subsequent RAM-backed Linux validation

The exact `a2b3b44a74f1b85c66ecaf93292c5a5e07707ff4` source archive
(SHA-256 `eea1adfc609352e0745a7e073b822e89f3c83f74f84433a7c142078620bdd540`)
was streamed into task-owned Linux aarch64 containers with `/work` on tmpfs
and no host bind mounts. Cargo, Rustup, pip, the Python environment and build
outputs stayed in RAM; reports and the wheel were streamed back to T7.
Both task containers were then removed.

With the pinned tool versions above, 21 Rust unit tests and six integration
tests passed, one non-asserting diagnostic remained ignored, and strict Clippy
passed. The release Linux wheel built successfully and all three installed-wheel
tests passed on Python 3.15.1 from outside the source directory. This validates
the local Linux ARM64 scope; it is not a GitLab runner result.

Actual MR pipeline `2986002132` created both jobs, but each failed with
`ci_quota_exceeded` and `started_at=null`. Neither job executed.

Independent read-only review of exported snapshot `ci-a2b3b44` completed through
OpenCode `openrouter/poolside/laguna-s-2.1:free` (reasoning variant `low`) with
scoped approval and no concrete blockers. Task:
`cea-ci-free-review-r1-snapshot-20261025`. The CI YAML is unchanged after review;
this section records subsequent evidence only. Numerical source is unchanged
from the independently reviewed performance snapshot `4734f3d`.
