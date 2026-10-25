"""Exercise an installed wheel from outside its source tree."""

from __future__ import annotations

import argparse
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--installed", type=Path, required=True)
args = parser.parse_args()
sys.path.insert(0, str(args.installed.resolve()))
import cea_rs


class InstalledWheelTests(unittest.TestCase):
    """Check public failures, invariants and repeated mixed workloads."""

    def test_install_and_bad_inputs(self) -> None:
        """Verify the import path and ensure failed solves expose no state."""
        self.assertTrue(Path(cea_rs.__file__).is_relative_to(args.installed.resolve()))
        solver = cea_rs.Solver(["H2", "O2"])
        for weights in ([], [0, 0], [-1, 1], [float("nan"), 1], [float("inf"), 1]):
            with self.assertRaises(ValueError):
                solver.solve("TP", 3000, 1e6, weights)
        with self.assertRaises(cea_rs.DomainError):
            cea_rs.Reactant("N2O").thermo(298.15)
        with self.assertRaises(cea_rs.DomainError):
            cea_rs.Reactant("RP-1").thermo(500)
        with self.assertRaises(ValueError):
            cea_rs.Solver(["H2"], only=[])
        with self.assertRaises(NotImplementedError):
            solver.solve("shock", 3000, 1e6, [1, 8])
        for options in ({"iac": False}, {"transport": True}):
            with self.assertRaises(NotImplementedError):
                solver.rocket([1, 6], [298.15] * 2, 3e6, 10, **options)
        impossible = cea_rs.Solver(["H2", "O2"], only=["H2O"])
        with self.assertRaises(cea_rs.ConvergenceError):
            impossible.solve("TP", 3000, 1e5, [1, 1])

    def test_interleaving_and_threads(self) -> None:
        """Keep cold/hot and frozen/equilibrium calls deterministic and isolated."""
        solver = cea_rs.Solver(["H2", "O2"])
        other = cea_rs.Solver(["N2O4(L)", "N2H4(L)", "C2H8N2(L),UDMH"])

        def hot(flow: str) -> tuple[float, ...]:
            """Return independent chamber, throat, exit and performance values."""
            r = solver.rocket([1, 6], [298.15] * 2, 3e6, 10, flow=flow)
            return (
                r.chamber.state.temperature_k,
                r.throat.state.temperature_k,
                r.exit.state.temperature_k,
                r.cstar_m_s,
                r.vacuum_velocity_m_s,
            )

        flows = ["equilibrium", "frozen_at_chamber", "frozen_at_throat"]
        baseline = {flow: hot(flow) for flow in flows}
        for flow in reversed(flows):
            solver.solve("TP", 500, 1e5, [1, 8])
            other.rocket(
                [1.6, 0.5, 0.5], [298.15] * 3, 848000, 45.6, flow="frozen_at_throat"
            )
            self.assertEqual(hot(flow), baseline[flow])
        with ThreadPoolExecutor(max_workers=3) as pool:
            for flow, values in zip(flows, pool.map(hot, flows), strict=True):
                self.assertEqual(values, baseline[flow])

    def test_nozzle_identities(self) -> None:
        """Check sonic, entropy, energy, mass and explicit frozen diagnostics."""
        solver = cea_rs.Solver(["H2", "O2"])
        for flow in ("equilibrium", "frozen_at_chamber", "frozen_at_throat"):
            r = solver.rocket([1, 6], [298.15] * 2, 3e6, 10, flow=flow)
            self.assertIsNone(r.chamber.area_ratio)
            self.assertAlmostEqual(r.throat.mach, 1.0, places=5)
            g = r.throat.state.density_kg_m3 * r.throat.velocity_m_s
            self.assertAlmostEqual(
                g / (r.exit.state.density_kg_m3 * r.exit.velocity_m_s * 10), 1, places=5
            )
            self.assertAlmostEqual(r.cstar_m_s, 3e6 / g, places=7)
            self.assertAlmostEqual(
                r.vacuum_velocity_m_s,
                r.exit.velocity_m_s + r.exit.state.pressure_pa * 10 / g,
                places=7,
            )
            if flow != "equilibrium":
                self.assertIsNone(r.exit.state.cp_equilibrium_j_kg_k)
            self.assertAlmostEqual(
                sum(r.exit.state.mass_fractions.values()), 1, places=10
            )
            self.assertAlmostEqual(
                r.chamber.state.enthalpy_j_kg,
                r.exit.state.enthalpy_j_kg + r.exit.velocity_m_s**2 / 2,
                places=5,
            )


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=2)
