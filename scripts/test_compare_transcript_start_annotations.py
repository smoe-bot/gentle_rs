#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import compare_transcript_start_annotations as target


class CompareTranscriptStartsTests(unittest.TestCase):
    def fixture(self, root: Path, secondary_rows=None):
        primary = root / "primary.json"
        catalog = root / "catalog.json"
        annotation = root / "annotation.gff"
        transcripts = root / "transcripts.json"
        manifest = root / "manifest.json"
        cli = root / "gentle_cli"
        cli.write_bytes(b"synthetic executable")
        annotation.write_text("##gff-version 3\n", encoding="utf-8")
        rows = secondary_rows or [
            {"chromosome":"NC_000001.11","transcript_id":"NM_PLUS.1","gene_id":"1",
             "gene_name":"PLUS","strand":"+","transcript_start_1based":1005,
             "transcript_end_1based":1500,"exons_1based":[[1005,1500]],"cds_1based":[]},
            {"chromosome":"NC_000002.12","transcript_id":"NM_MINUS.1","gene_id":"2",
             "gene_name":"MINUS","strand":"-","transcript_start_1based":1500,
             "transcript_end_1based":2010,"exons_1based":[[1500,2010]],"cds_1based":[]},
        ]
        transcripts.write_text(json.dumps(rows), encoding="utf-8")
        manifest.write_text(json.dumps({
            "genome_id":"Human GRCh38 NCBI RefSeq GCF_000001405.40",
            "annotation_path":str(annotation),
            "annotation_source":"https://example.invalid/refseq.gff.gz",
            "annotation_sha1":hashlib.sha1(annotation.read_bytes()).hexdigest(),
            "transcript_index_path":str(transcripts),
        }), encoding="utf-8")
        catalog.write_text(json.dumps({
            "Human GRCh38 Ensembl 116": {
                "description":"primary", "ensembl_template":{"file_stem":"Homo_sapiens.GRCh38","release":116}},
            "Human GRCh38 NCBI RefSeq GCF_000001405.40": {
                "description":"secondary", "ncbi_assembly_accession":"GCF_000001405.40",
                "ncbi_assembly_name":"GRCh38.p14"},
        }), encoding="utf-8")
        primary.write_text(json.dumps({
            "schema":target.TSS_REPORT_SCHEMA,
            "reference":{"genome_id":"Human GRCh38 Ensembl 116","assembly":"GRCh38"},
            "windows":[
                {"selected":True,"record":{"gene_symbol":"PLUS","promoter_id":"p1",
                 "geometry":{"chromosome":"1","strand":"+","tss_1based":1000},"transcripts":["ENSTP1"]}},
                {"selected":True,"record":{"gene_symbol":"MINUS","promoter_id":"m1",
                 "geometry":{"chromosome":"2","strand":"-","tss_1based":2000},"transcripts":["ENSTM1"]}},
            ],
        }), encoding="utf-8")
        return argparse.Namespace(
            primary_tss_report=primary, secondary_manifest=manifest, catalog=catalog,
            gene=["PLUS","MINUS"], output_report=root/"out.json", output_dir=root/"pages",
            output_receipt=root/"receipt.json", producer_revision="test-revision", gentle_cli=cli)

    def test_plus_and_minus_deltas_are_transcript_oriented(self):
        with tempfile.TemporaryDirectory() as raw:
            args = self.fixture(Path(raw))
            report = target.build_report(args)
            plus, minus = report["genes"]
            self.assertEqual(plus["primary_clusters"][0]["nearest_secondary"]["delta_bp_transcript_direction"], 5)
            self.assertEqual(minus["primary_clusters"][0]["nearest_secondary"]["delta_bp_transcript_direction"], -10)
            self.assertEqual(plus["secondary_transcript_count"], 1)
            self.assertEqual(minus["secondary_transcript_count"], 1)

    def test_outputs_bind_report_and_aligned_pages(self):
        with tempfile.TemporaryDirectory() as raw:
            args = self.fixture(Path(raw))
            report = target.build_report(args)
            target.write_outputs(args, report)
            receipt = json.loads(args.output_receipt.read_text())
            self.assertEqual(receipt["schema"], target.RECEIPT_SCHEMA)
            self.assertEqual(len(receipt["pages"]), 2)
            svg = (args.output_dir / "PLUS_Ensembl_RefSeq_TSS_comparison.svg").read_text()
            self.assertIn('data-gentle-plot-left="255"', svg)
            self.assertIn("+5 bp", svg)

    def test_annotation_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as raw:
            args = self.fixture(Path(raw))
            manifest = json.loads(args.secondary_manifest.read_text())
            manifest["annotation_sha1"] = "0" * 40
            args.secondary_manifest.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "SHA-1"):
                target.build_report(args)

    def test_assembly_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as raw:
            args = self.fixture(Path(raw))
            catalog = json.loads(args.catalog.read_text())
            catalog["Human GRCh38 NCBI RefSeq GCF_000001405.40"]["ncbi_assembly_name"] = "T2T-CHM13v2.0"
            args.catalog.write_text(json.dumps(catalog))
            with self.assertRaisesRegex(ValueError, "same assembly"):
                target.build_report(args)


if __name__ == "__main__":
    unittest.main()
