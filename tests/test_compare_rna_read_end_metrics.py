import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compare_rna_read_end_metrics.py"


class CompareRnaReadEndMetricsTests(unittest.TestCase):
    def write_fixture(self, path: Path) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "metric", "value_bp", "observation_status"])
            for metric_index, metric in enumerate(("read_length_bp", "utr5_reach_bp", "utr3_reach_bp")):
                for sample_index in range(6):
                    for value in (100, 200, 300, 400):
                        writer.writerow([f"barcode{13 + sample_index}", metric, value + 10 * metric_index + sample_index, "observed"])
                    writer.writerow([f"barcode{13 + sample_index}", metric, "", "missing"])

    def test_six_samples_three_metrics_produce_all_fifteen_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "metrics.tsv", root / "out"
            self.write_fixture(source)
            subprocess.run([sys.executable, str(SCRIPT), str(source), str(output)], check=True)
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["schema"], "gentle.rna_read_metric_comparison.v1")
            self.assertEqual(report["pair_count_per_metric"], 15)
            self.assertEqual(len(report["pairs"]), 45)
            with (output / "qq_pairs.tsv").open(encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in handle) - 1, 15 * 3 * 101)
            summary = (output / "summary.tsv").read_text(encoding="utf-8")
            self.assertIn("barcode13\tread_length_bp\t4\t0\t0\t1\t", summary)
            for metric in ("read_length_bp", "utr5_reach_bp", "utr3_reach_bp"):
                self.assertIn("gray diagonal is identity", (output / f"qq_{metric}.svg").read_text(encoding="utf-8"))

    def test_missing_or_censored_value_cannot_masquerade_as_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "bad.tsv"
            source.write_text(
                "sample_id\tmetric\tvalue_bp\tobservation_status\n"
                "a\tread_length_bp\t100\tobserved\n"
                "a\tread_length_bp\t0\tmissing\n",
                encoding="utf-8",
            )
            completed = subprocess.run([sys.executable, str(SCRIPT), str(source), str(root / "out")], text=True, capture_output=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("never encoded as zero", completed.stderr)


if __name__ == "__main__":
    unittest.main()
