"""Compare two installed versions in alternating separate-process rounds."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path


def main() -> None:
    """Retain raw same-input timings and pair ratios, avoiding stale baselines."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    work = args.output.resolve().parent / "paired-rounds"
    work.mkdir(parents=True, exist_ok=True)
    rounds = []
    for round_number in range(3):
        results = {}
        order = ["baseline", "candidate"]
        if round_number % 2:
            order.reverse()
        for label in order:
            destination = work / f"{round_number}-{label}.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).with_name("benchmark.py")),
                    "--installed",
                    str(getattr(args, label).resolve()),
                    "--output",
                    str(destination),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            (work / f"{round_number}-{label}.log").write_text(
                completed.stdout + completed.stderr
            )
            results[label] = json.loads(destination.read_text())
        rounds.append(results)
    cases = []
    for i, baseline in enumerate(rounds[0]["baseline"]["cases"]):
        ratios = []
        for results in rounds:
            old = results["baseline"]["cases"][i]
            new = results["candidate"]["cases"][i]
            assert old["reactants"] == new["reactants"] and old["mode"] == new["mode"]
            ratios.append(
                old["warm"]["rust"]["median_ms"] / new["warm"]["rust"]["median_ms"]
            )
        cases.append(
            {
                "reactants": baseline["reactants"],
                "mode": baseline["mode"],
                "paired_baseline_over_candidate": ratios,
                "median_speedup": statistics.median(ratios),
            }
        )
    report = {
        "rounds": rounds,
        "cases": cases,
        "conditions": "Three alternating baseline/candidate process rounds, same Python interpreter and benchmark inputs. Each has five warmups/50 warm samples plus five fresh-process cold samples per backend/mode. OS caches not flushed; CPU scheduling/thermal conditions uncontrolled. TP is unchanged and serves as a noise check.",
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for case in cases:
        print(
            case["mode"],
            case["reactants"],
            round(case["median_speedup"], 3),
            flush=True,
        )


if __name__ == "__main__":
    main()
