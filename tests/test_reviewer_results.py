import csv
import json
import statistics
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "runs" / "reviewer"


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class ReviewerResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = read_csv("grid_and_filter_sensitivity.csv")
        cls.packing = read_csv("packing_sensitivity.csv")

    def test_result_counts_match_summary(self):
        summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary, {"grid_rows": 42, "packing_rows": 140})
        self.assertEqual(len(self.grid), summary["grid_rows"])
        self.assertEqual(len(self.packing), summary["packing_rows"])

    def test_grid_convergence_claims(self):
        rows = {
            (row["framework"], int(row["grid_n"]), int(row["n_min"])): row
            for row in self.grid
        }

        def accessible(framework, grid_n):
            return float(rows[framework, grid_n, 4]["accessible_cell_fraction"])

        fau_pp = abs(accessible("FAU", 48) - accessible("FAU", 96)) * 100
        fau_relative_percent = fau_pp / (accessible("FAU", 96) * 100) * 100
        self.assertEqual(round(fau_pp, 2), 0.45)
        self.assertEqual(round(fau_relative_percent, 2), 1.46)

        others = ["LTL", "ERI", "MOF-5", "UiO-66", "COF-5"]
        max_other_pp = max(
            abs(accessible(name, 48) - accessible(name, 96)) * 100
            for name in others
        )
        self.assertEqual(round(max_other_pp, 2), 0.39)
        self.assertEqual(
            [int(rows["FAU", n, 4]["blocked_clusters"]) for n in (32, 48, 64)],
            [8, 8, 8],
        )

    def test_default_iteration_diagnostics(self):
        expected_means = {
            ("LTL", "EC"): 8.8,
            ("UiO-66", "EC"): 3.2,
            ("UiO-66", "LiPF6"): 1.4,
            ("COF-5", "EC"): 3.4,
        }
        for case, expected in expected_means.items():
            rows = [
                row
                for row in self.packing
                if (row["framework"], row["guest"]) == case
                and row["cutoff_A"] == "2.2"
                and row["blocked_penalty"] == "100.0"
                and row["max_iters"] == "12000"
            ]
            self.assertEqual(len(rows), 5)
            self.assertTrue(all(int(row["inaccessible"]) == 0 for row in rows))
            self.assertEqual(
                statistics.mean(int(row["total_violations"]) for row in rows),
                expected,
            )

    def test_cutoff_and_penalty_claims(self):
        fau_co2_cutoffs = [
            row
            for row in self.packing
            if row["framework"] == "FAU"
            and row["guest"] == "CO2"
            and row["blocked_penalty"] == "100.0"
            and row["max_iters"] == "30000"
        ]
        self.assertTrue(fau_co2_cutoffs)
        self.assertTrue(all(int(row["total_violations"]) == 0 for row in fau_co2_cutoffs))

        h2o_26 = [
            row
            for row in self.packing
            if row["framework"] == "FAU"
            and row["guest"] == "H2O"
            and row["cutoff_A"] == "2.6"
        ]
        self.assertEqual(len(h2o_26), 3)
        self.assertEqual(sum(int(row["total_violations"]) > 0 for row in h2o_26), 1)

        penalty_rows = [
            row
            for row in self.packing
            if row["framework"] == "FAU"
            and row["guest"] == "CO2"
            and row["cutoff_A"] == "2.2"
            and row["max_iters"] == "30000"
        ]
        tested_penalties = {float(row["blocked_penalty"]) for row in penalty_rows}
        self.assertTrue({0.1, 1.0, 10.0, 100.0, 1000.0}.issubset(tested_penalties))
        self.assertTrue(all(int(row["total_violations"]) == 0 for row in penalty_rows))


if __name__ == "__main__":
    unittest.main()
