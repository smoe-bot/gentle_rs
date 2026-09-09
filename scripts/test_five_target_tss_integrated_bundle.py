#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "docs/examples/regulatory_region_comparison/five_target_tss_integrated"
FASTAS = ROOT / "docs/examples/regulatory_region_comparison/target_tss_fastas/bundle"
GENES = {"CD44", "TGFB1", "SERPINE1", "PATZ1", "TP73"}


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


class FiveTargetTssIntegratedBundleTests(unittest.TestCase):
    def test_checksum_inventory_is_complete(self) -> None:
        rows = [line.split("  ", 1) for line in (BUNDLE / "SHA256SUMS").read_text().splitlines()]
        recorded = {name: value for value, name in rows}
        actual = {path.name for path in BUNDLE.iterdir() if path.is_file() and path.name != "SHA256SUMS"}
        self.assertEqual(set(recorded), actual)
        self.assertTrue(all(digest(BUNDLE / name) == value for name, value in recorded.items()))

    def test_reports_have_complete_track_and_layout_contract(self) -> None:
        for gene in GENES:
            stem = f"{gene}_luciferase_planning_EnsemblReg_30TF_with_TSS_similarity"
            svg = (BUNDLE / f"{stem}.svg").read_text()
            self.assertEqual(svg.count("JASPAR PWM"), 30)
            self.assertIn("TFAP2C JASPAR PWM (MA0524.3)", svg)
            self.assertIn("TFAP2C JASPAR PWM (MA0814.3)", svg)
            self.assertIn("TFAP2C JASPAR PWM (MA0815.1)", svg)
            self.assertIn('data-gentle-tss-stretch-backgrounds="true"', svg)
            self.assertIn('data-gentle-tss-stretch-band=', svg)
            self.assertRegex(svg, r"chr(?:1|7|11|19|22):[0-9]+")
            self.assertIn("TSS-local Ensembl-feature promoterome similarity", svg)
            self.assertTrue((BUNDLE / f"{stem}.pdf").read_bytes().startswith(b"%PDF"))
            self.assertTrue((BUNDLE / f"{stem}.png").read_bytes().startswith(b"\x89PNG"))
            receipt = json.loads((BUNDLE / f"{stem}.receipt.json").read_text())
            for filename, expected in receipt["outputs"].items():
                self.assertEqual("sha256:" + digest(BUNDLE / filename), expected)

    def test_scientific_inputs_and_fastas_cover_five_targets(self) -> None:
        candidates = json.loads((BUNDLE / "candidate_regions.json").read_text())
        self.assertEqual({row["gene_query"] for row in candidates["regions"]}, GENES)
        self.assertEqual(len(candidates["regions"]), 26)
        lanes = json.loads((BUNDLE / "cutrun_lane_validation.json").read_text())
        self.assertEqual({row["gene"] for row in lanes["genes"]}, GENES)
        self.assertEqual(lanes["total_lane_count"], 60)
        self.assertEqual(lanes["total_rendered_interval_count"], 35211)
        panel = json.loads((BUNDLE / "jaspar_30_track_panel.json").read_text())["tracks"]
        self.assertEqual(len(panel), 30)
        self.assertEqual(len({track["factors"][0]["factor_id"] for track in panel}), 28)
        for gene in GENES:
            fasta = FASTAS / f"{gene}_TSS_minus500_plus200.fasta"
            text = fasta.read_text()
            self.assertGreater(text.count(">"), 0)
            lengths = []
            for record in re.split(r"(?m)^>", text)[1:]:
                lines = record.splitlines()
                lengths.append(len("".join(lines[1:])))
            self.assertEqual(set(lengths), {701})


if __name__ == "__main__":
    unittest.main()
