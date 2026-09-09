#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
import unittest

from scripts import build_prepared_reference_gene_entry as builder


class PreparedReferenceGeneEntryTests(unittest.TestCase):
    def fixture(self, root: Path, strand: str) -> argparse.Namespace:
        fasta = root / "reference.fa"
        sequence = "ACGT" * 100
        fasta.write_text(">1\n" + "\n".join(sequence[i:i + 20] for i in range(0, len(sequence), 20)) + "\n")
        fai = root / "reference.fa.fai"
        fai.write_text("1\t400\t3\t20\t21\n")
        transcripts = root / "transcripts.json"
        transcripts.write_text(json.dumps([{
            "chromosome": "1", "transcript_id": "ENST1", "gene_id": "ENSG1",
            "gene_name": "GENE", "strand": strand, "transcript_start_1based": 101,
            "transcript_end_1based": 160, "exons_1based": [[101, 120], [141, 160]],
            "cds_1based": [[105, 120], [141, 154]],
        }]))
        return argparse.Namespace(transcripts=transcripts, fasta=fasta, fai=fai, gene="GENE",
            species="homo_sapiens", assembly="GRCh38", annotation_release="Ensembl 116",
            flank_5prime_bp=10, flank_3prime_bp=5, entry_id="GENE_entry", transcript_id=None)

    def test_plus_and_minus_are_transcript_oriented(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            plus = builder.build(self.fixture(root, "+"))
            minus = builder.build(self.fixture(root, "-"))
            self.assertEqual((plus["sequence_genomic_start_1based"], plus["sequence_genomic_end_1based"]), (91, 165))
            self.assertEqual((minus["sequence_genomic_start_1based"], minus["sequence_genomic_end_1based"]), (96, 170))
            genomic, _ = builder.fetch_region(root / "reference.fa", root / "reference.fa.fai", "1", 96, 170)
            self.assertEqual(minus["sequence"], genomic.translate(builder.DNA_COMPLEMENT)[::-1])
            self.assertEqual(plus["transcripts"][0]["translation"]["length_aa"], 10)

    def test_rejects_mixed_gene_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.fixture(root, "+")
            rows = json.loads(args.transcripts.read_text())
            rows.append({**rows[0], "transcript_id": "ENST2", "chromosome": "2"})
            args.transcripts.write_text(json.dumps(rows))
            with self.assertRaisesRegex(RuntimeError, "disagree"):
                builder.build(args)


if __name__ == "__main__":
    unittest.main()
