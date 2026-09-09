#!/usr/bin/env python3
"""Build a deterministic offline Ensembl-gene entry from a prepared GENtle reference."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

try:
    from .prepare_regulatory_region_indexes import require
except ImportError:
    from prepare_regulatory_region_indexes import require


DNA_COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_fai(path: Path, chromosome: str) -> tuple[int, int, int, int]:
    with path.open(encoding="ascii") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == chromosome:
                return tuple(map(int, fields[1:5]))  # type: ignore[return-value]
    raise RuntimeError(f"chromosome absent from FASTA index: {chromosome}")


def fetch_region(fasta: Path, fai: Path, chromosome: str, start_1based: int,
                 end_1based: int) -> tuple[str, int]:
    length, offset, line_bases, line_width = read_fai(fai, chromosome)
    require(1 <= start_1based <= end_1based <= length, "requested sequence region is out of bounds")
    start0 = start_1based - 1
    end0 = end_1based
    byte_start = offset + (start0 // line_bases) * line_width + start0 % line_bases
    byte_end = offset + (end0 // line_bases) * line_width + end0 % line_bases
    with fasta.open("rb") as handle:
        handle.seek(byte_start)
        payload = handle.read(byte_end - byte_start + line_width)
    sequence = b"".join(payload.splitlines()).decode("ascii")[:end0 - start0].upper()
    require(len(sequence) == end0 - start0 and set(sequence) <= set("ACGTRYMKSWBDHVN"),
            "prepared-reference sequence has invalid length or alphabet")
    return sequence, length


def build(args: argparse.Namespace) -> dict:
    transcripts_path = args.transcripts.resolve(strict=True)
    fasta = args.fasta.resolve(strict=True)
    fai = args.fai.resolve(strict=True)
    records = json.loads(transcripts_path.read_text(encoding="utf-8"))
    require(isinstance(records, list), "transcript index must be a JSON array")
    gene_records = [row for row in records if row.get("gene_name") == args.gene]
    requested = set(args.transcript_id or [])
    selected = ([row for row in gene_records if row.get("transcript_id") in requested]
                if requested else gene_records)
    require(selected, f"gene absent from transcript index: {args.gene}")
    require(not requested or {row["transcript_id"] for row in selected} == requested,
            "one or more requested transcripts are absent from the target gene")
    gene_ids = {row.get("gene_id") for row in selected}
    chromosomes = {row.get("chromosome") for row in selected}
    strands = {row.get("strand") for row in selected}
    require(len(gene_ids) == len(chromosomes) == len(strands) == 1,
            "gene records disagree on gene ID, chromosome, or strand")
    gene_id = next(iter(gene_ids))
    chromosome = next(iter(chromosomes))
    strand_token = next(iter(strands))
    require(isinstance(gene_id, str) and gene_id and isinstance(chromosome, str) and chromosome,
            "gene identity is incomplete")
    require(strand_token in {"+", "-"}, "gene strand must be '+' or '-'")
    strand = 1 if strand_token == "+" else -1
    gene_start = min(int(row["transcript_start_1based"]) for row in selected)
    gene_end = max(int(row["transcript_end_1based"]) for row in selected)
    contig_length = read_fai(fai, chromosome)[0]
    if strand == 1:
        sequence_start = max(1, gene_start - args.flank_5prime_bp)
        sequence_end = min(contig_length, gene_end + args.flank_3prime_bp)
    else:
        sequence_start = max(1, gene_start - args.flank_3prime_bp)
        sequence_end = min(contig_length, gene_end + args.flank_5prime_bp)
    sequence, _ = fetch_region(fasta, fai, chromosome, sequence_start, sequence_end)
    if strand == -1:
        sequence = sequence.translate(DNA_COMPLEMENT)[::-1]

    transcript_summaries = []
    for row in sorted(selected, key=lambda item: item["transcript_id"]):
        exons = [{
            "exon_id": f"{row['transcript_id']}_exon_{index}", "exon_version": None,
            "start_1based": int(span[0]), "end_1based": int(span[1]),
            "strand": strand, "seq_region_name": chromosome,
        } for index, span in enumerate(row.get("exons_1based", []), 1)]
        cds = row.get("cds_1based", [])
        translation = None
        if cds:
            coding_nt = sum(int(span[1]) - int(span[0]) + 1 for span in cds)
            translation = {
                "translation_id": f"{row['transcript_id']}_translation", "translation_version": None,
                "length_aa": coding_nt // 3, "genomic_start_1based": min(int(x[0]) for x in cds),
                "genomic_end_1based": max(int(x[1]) for x in cds),
            }
        transcript_summaries.append({
            "transcript_id": row["transcript_id"], "transcript_version": None,
            "display_name": None, "biotype": "protein_coding" if cds else "noncoding_or_no_cds",
            "start_1based": int(row["transcript_start_1based"]),
            "end_1based": int(row["transcript_end_1based"]), "strand": strand,
            "is_canonical": None, "gencode_primary": None, "translation": translation, "exons": exons,
        })

    binding = {
        "transcript_index_sha256": sha256_file(transcripts_path),
        "sequence_fai_sha256": sha256_file(fai), "assembly": args.assembly,
        "annotation_release": args.annotation_release,
    }
    return {
        "schema": "gentle.ensembl_gene_entry.v1", "entry_id": args.entry_id,
        "gene_id": gene_id, "gene_version": None, "gene_symbol": args.gene,
        "gene_display_name": args.gene, "species": args.species,
        "assembly_name": args.assembly, "biotype": "gene", "strand": strand,
        "seq_region_name": chromosome, "genomic_start_1based": gene_start,
        "genomic_end_1based": gene_end, "sequence_genomic_start_1based": sequence_start,
        "sequence_genomic_end_1based": sequence_end, "flank_5prime_bp": args.flank_5prime_bp,
        "flank_3prime_bp": args.flank_3prime_bp, "sequence": sequence,
        "sequence_length": len(sequence), "transcripts": transcript_summaries,
        "aliases": [], "source": "GENtle prepared Human GRCh38 Ensembl 116 reference",
        "source_query": args.gene, "imported_at_unix_ms": 0,
        "lookup_source_url": "", "sequence_source_url": "",
        "raw_lookup_json": json.dumps(binding, sort_keys=True, separators=(",", ":")),
        "raw_sequence_json": json.dumps({**binding, "chromosome": chromosome,
            "start_1based": sequence_start, "end_1based": sequence_end,
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest()},
            sort_keys=True, separators=(",", ":")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transcripts", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--transcript-id", action="append")
    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--annotation-release", default="Ensembl 116")
    parser.add_argument("--flank-5prime-bp", type=int, default=5000)
    parser.add_argument("--flank-3prime-bp", type=int, default=1000)
    parser.add_argument("--entry-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.gene and args.gene.upper() == args.gene, "--gene must be an uppercase symbol")
    require(args.flank_5prime_bp >= 0 and args.flank_3prime_bp >= 0, "flanks must be non-negative")
    payload = build(args)
    require(not args.output.exists(), "output already exists")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "gene": args.gene,
                      "transcripts": len(payload["transcripts"]), "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
