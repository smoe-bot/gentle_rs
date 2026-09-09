#!/usr/bin/env python3
"""Focused tests for target-gene TSS FASTA export."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts import export_target_tss_fastas as exporter
from scripts.prepare_transcript_promoterome import RECEIPT_SCHEMA, window_id


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


class ExportTargetTssFastasTests(unittest.TestCase):
    def make_promoterome(self, root: Path, *, duplicate_gene_id: bool = False) -> Path:
        bundle = root / "promoterome"
        bundle.mkdir()
        rows = []
        mappings = []
        fasta = []
        specs = [
            ("1", 10000, "+", "A" * 2201, ["ENST1", "ENST2"]),
            ("2", 20000, "-", "C" * 1500 + "G" * 701, ["ENST3"]),
        ]
        for chromosome, tss, strand, sequence, transcripts in specs:
            start = tss - 2000 - 1 if strand == "+" else tss - 200 - 1
            end = tss + 200 if strand == "+" else tss + 2000
            promoter = window_id(chromosome, start, end, strand)
            rows.append([promoter, chromosome, start, end, strand, tss, "false", 1, len(transcripts)])
            fasta.append(f">{promoter}\n{sequence}\n")
            for index, transcript in enumerate(transcripts):
                gene_id = "ENSG2" if duplicate_gene_id and index else "ENSG1"
                mappings.append([promoter, gene_id, "GENE1", transcript])
        with (bundle / "promoter_windows.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["promoter_id", "chromosome", "start_0based", "end_0based_exclusive",
                             "strand", "tss_1based", "boundary_clipped", "gene_count", "transcript_count"])
            writer.writerows(rows)
        with (bundle / "promoter_transcripts.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["promoter_id", "gene_id", "gene_name", "transcript_id"])
            writer.writerows(mappings)
        (bundle / "promoter_windows.fa").write_text("".join(fasta), encoding="ascii")
        artifacts = {name: "sha256:" + digest(bundle / name) for name in
                     ("promoter_windows.tsv", "promoter_transcripts.tsv", "promoter_windows.fa")}
        (bundle / "receipt.json").write_text(json.dumps({
            "schema": RECEIPT_SCHEMA, "dataset_id": "synthetic", "genome_id": "Human GRCh38 Ensembl 116",
            "upstream_bp": 2000, "downstream_bp": 200, "unique_promoter_window_count": 2,
            "included_transcript_count": 3, "sequence_orientation": "transcript_5prime_to_3prime_via_bedtools_strand",
            "artifacts": artifacts,
        }), encoding="utf-8")
        return bundle

    def args(self, root: Path, promoterome: Path) -> argparse.Namespace:
        return argparse.Namespace(promoterome=promoterome, gene=["GENE1"], upstream_bp=500,
                                  downstream_bp=200, expected_genome_id="Human GRCh38 Ensembl 116",
                                  assembly_id="GRCh38", source_revision="a" * 40, output=root / "out")

    @patch.object(subprocess, "check_output", return_value="a" * 40 + "\n")
    def test_exports_distinct_tss_records_and_shared_transcripts(self, _revision) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = exporter.export(self.args(root, self.make_promoterome(root)))
            self.assertEqual(manifest["files"][0]["record_count"], 2)
            text = (root / "out/GENE1_TSS_minus500_plus200.fasta").read_text()
            self.assertEqual(text.count(">"), 2)
            self.assertIn("transcripts=ENST1,ENST2", text)
            records = text.split(">")
            self.assertIn("\n" + "A" * 80, records[1])
            self.assertIn("\n" + "G" * 80, records[2])
            self.assertIn("strand=-", records[2])
            self.assertTrue((root / "out/SHA256SUMS").is_file())

    @patch.object(subprocess, "check_output", return_value="a" * 40 + "\n")
    def test_rejects_ambiguous_gene_symbol(self, _revision) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with self.assertRaisesRegex(RuntimeError, "exactly one"):
                exporter.export(self.args(root, self.make_promoterome(root, duplicate_gene_id=True)))

    @patch.object(subprocess, "check_output", return_value="a" * 40 + "\n")
    def test_rejects_wrong_genome_identity(self, _revision) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.args(root, self.make_promoterome(root))
            args.expected_genome_id = "Human GRCh37 Ensembl 116"
            with self.assertRaisesRegex(RuntimeError, "exactly match"):
                exporter.export(args)


if __name__ == "__main__":
    unittest.main()
