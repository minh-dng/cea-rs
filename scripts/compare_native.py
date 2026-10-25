"""Run native NASA comparison fixtures in isolated worker processes."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def main() -> None:
    """Isolate native crashes, retain each fixture result, and merge evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    outputs = args.output.resolve().parent / "native-cases"
    outputs.mkdir(parents=True, exist_ok=True)
    sections = [
        "equilibrium-0",
        "equilibrium-1",
        "condensed",
        "ions",
        *[f"rocket-{i}" for i in range(5)],
    ]
    reports = []
    failures = []
    for section in sections:
        output = outputs / f"{section}.json"
        command = [
            sys.executable,
            str(Path(__file__).with_name("_compare_case.py")),
            "--installed",
            str(args.installed.resolve()),
            "--output",
            str(output),
            "--section",
            section,
        ]
        try:
            result = subprocess.run(
                command, capture_output=True, text=True, timeout=180, check=False
            )
            (outputs / f"{section}.log").write_text(result.stdout + result.stderr)
            print(section, result.returncode, flush=True)
            if result.returncode != 0:
                failures.append({"section": section, "returncode": result.returncode})
            else:
                reports.append(json.loads(output.read_text()))
        except subprocess.TimeoutExpired:
            failures.append({"section": section, "error": "timeout after 180 s"})
    metadata = (
        {k: v for k, v in reports[0].items() if k not in ("cases", "all_passed")}
        if reports
        else {}
    )
    merged = {
        **metadata,
        "process_isolated": True,
        "cases": [c for r in reports for c in r["cases"]],
        "worker_failures": failures,
        "all_passed": not failures
        and len(reports) == len(sections)
        and all(r["all_passed"] for r in reports),
    }
    args.output.write_text(json.dumps(merged, indent=2, allow_nan=False) + "\n")
    if not merged["all_passed"]:
        raise SystemExit("Native comparisons failed; see isolated worker logs")


if __name__ == "__main__":
    main()
