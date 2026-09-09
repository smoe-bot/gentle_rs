#!/usr/bin/env python3
"""Export one receipt-bound, transcript-oriented TSS FASTA per target gene."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

try:
    from .prepare_transcript_promoterome import validate_promoterome, window_id
    from .prepare_regulatory_region_indexes import require
except ImportError:
    from prepare_transcript_promoterome import validate_promoterome, window_id
    from prepare_regulatory_region_indexes import require


SCHEMA = "gentle.target_tss_fasta_export.v1"
DNA = set("ACGTRYSWKMBDHVN")


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    found: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                current = name if name in wanted else None
                if current is not None:
                    require(current not in found, f"duplicate promoter FASTA record: {current}")
                    found[current] = []
            elif current is not None and line:
                found[current].append(line.upper())
    require(set(found) == wanted, f"missing promoter FASTA records: {sorted(wanted - set(found))}")
    return {name: "".join(parts) for name, parts in found.items()}


def chromosome_key(value: str) -> tuple[int, int | str]:
    token = value.removeprefix("chr")
    return (0, int(token)) if token.isdigit() else (1, token)


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def export(args: argparse.Namespace) -> dict[str, Any]:
    require(args.upstream_bp >= 0 and args.downstream_bp >= 0,
            "requested TSS offsets must be non-negative")
    genes = list(dict.fromkeys(args.gene))
    require(genes and all(gene and gene.upper() == gene for gene in genes),
            "--gene values must be non-empty uppercase gene symbols")
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "output directory must be absent or empty")
    promoterome = args.promoterome.resolve(strict=True)
    receipt = validate_promoterome(promoterome)
    require(receipt["genome_id"] == args.expected_genome_id,
            "promoterome genome_id does not exactly match --expected-genome-id")
    source_upstream = receipt["upstream_bp"]
    source_downstream = receipt["downstream_bp"]
    require(args.upstream_bp <= source_upstream and args.downstream_bp <= source_downstream,
            "requested window is not contained in the prepared promoterome window")

    repo = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    require(args.source_revision == revision,
            "--source-revision must be the current checkout's full HEAD SHA")

    window_rows = load_tsv(promoterome / "promoter_windows.tsv")
    mapping_rows = load_tsv(promoterome / "promoter_transcripts.tsv")
    require(len(window_rows) == receipt["unique_promoter_window_count"],
            "promoter window count disagrees with receipt")
    require(len(mapping_rows) == receipt["included_transcript_count"],
            "promoter transcript count disagrees with receipt")
    windows = {row["promoter_id"]: row for row in window_rows}
    require(len(windows) == len(window_rows), "duplicate promoter window identifiers")

    selected = [row for row in mapping_rows if row["gene_name"] in genes]
    by_gene: dict[str, list[dict[str, str]]] = {gene: [] for gene in genes}
    for row in selected:
        by_gene[row["gene_name"]].append(row)
    for gene, rows in by_gene.items():
        require(rows, f"target gene absent from promoterome: {gene}")
        gene_ids = {row["gene_id"] for row in rows}
        require(len(gene_ids) == 1 and "" not in gene_ids,
                f"target symbol does not resolve to exactly one non-empty gene ID: {gene}")
        require(len({row["transcript_id"] for row in rows}) == len(rows),
                f"duplicate transcript mapping for target gene: {gene}")
        require(all(row["promoter_id"] in windows for row in rows),
                f"target mapping refers to an absent promoter window: {gene}")

    wanted = {row["promoter_id"] for row in selected}
    source_sequences = extract_fasta(promoterome / "promoter_windows.fa", wanted)
    offset = source_upstream - args.upstream_bp
    length = args.upstream_bp + args.downstream_bp + 1
    output.mkdir(parents=True, exist_ok=True)
    manifest_files: list[dict[str, Any]] = []

    for gene in genes:
        mappings = by_gene[gene]
        grouped: dict[str, list[dict[str, str]]] = {}
        for row in mappings:
            grouped.setdefault(row["promoter_id"], []).append(row)
        records: list[dict[str, Any]] = []
        for promoter_id, members in grouped.items():
            window = windows[promoter_id]
            start0 = int(window["start_0based"])
            end0 = int(window["end_0based_exclusive"])
            tss = int(window["tss_1based"])
            strand = window["strand"]
            chromosome = window["chromosome"]
            require(strand in {"+", "-"} and str(window["boundary_clipped"]).lower() == "false",
                    f"selected promoter is clipped or has invalid strand: {promoter_id}")
            expected = ((tss - source_upstream - 1, tss + source_downstream) if strand == "+"
                        else (tss - source_downstream - 1, tss + source_upstream))
            require((start0, end0) == expected,
                    f"selected promoter geometry disagrees with its TSS: {promoter_id}")
            require(promoter_id == window_id(chromosome, start0, end0, strand),
                    f"selected promoter ID disagrees with its geometry: {promoter_id}")
            source_sequence = source_sequences[promoter_id]
            require(len(source_sequence) == end0 - start0 and set(source_sequence) <= DNA,
                    f"selected promoter sequence is invalid: {promoter_id}")
            sequence = source_sequence[offset:offset + length]
            require(len(sequence) == length, f"could not extract requested window: {promoter_id}")
            if strand == "+":
                genomic_start, genomic_end = tss - args.upstream_bp, tss + args.downstream_bp
            else:
                genomic_start, genomic_end = tss - args.downstream_bp, tss + args.upstream_bp
            transcripts = sorted({row["transcript_id"] for row in members})
            gene_id = members[0]["gene_id"]
            sequence_hash = sha256_text(sequence)
            header = (
                f"{gene}|gene_id={gene_id}|promoter_id={promoter_id}|assembly={args.assembly_id}"
                f"|chromosome={chromosome}|strand={strand}|tss_1based={tss}"
                f"|genomic_1based={genomic_start}-{genomic_end}"
                f"|window=minus{args.upstream_bp}_plus{args.downstream_bp}"
                f"|orientation=transcript_5prime_to_3prime|transcripts={','.join(transcripts)}"
                f"|sequence_sha256={sequence_hash}"
            )
            records.append({
                "promoter_id": promoter_id, "gene_id": gene_id, "chromosome": chromosome,
                "strand": strand, "tss_1based": tss, "genomic_start_1based": genomic_start,
                "genomic_end_1based": genomic_end, "transcript_ids": transcripts,
                "sequence_length_bp": len(sequence), "sequence_sha256": sequence_hash,
                "header": header, "sequence": sequence,
            })
        records.sort(key=lambda row: (chromosome_key(row["chromosome"]), row["tss_1based"], row["strand"]))
        filename = f"{gene}_TSS_minus{args.upstream_bp}_plus{args.downstream_bp}.fasta"
        fasta = output / filename
        fasta.write_text("".join(f">{row['header']}\n{wrap_fasta(row['sequence'])}\n" for row in records),
                         encoding="ascii")
        manifest_files.append({
            "gene_symbol": gene, "gene_id": records[0]["gene_id"], "filename": filename,
            "sha256": sha256_file(fasta), "record_count": len(records),
            "records": [{key: value for key, value in row.items() if key not in {"header", "sequence"}}
                        for row in records],
        })

    sums = output / "SHA256SUMS"
    sums.write_text("".join(f"{item['sha256']}  {item['filename']}\n" for item in manifest_files),
                    encoding="ascii")
    manifest = {
        "schema": SCHEMA, "source_revision": revision,
        "producer_sha256": sha256_file(Path(__file__)),
        "source": {
            "dataset_id": receipt["dataset_id"], "genome_id": receipt["genome_id"],
            "promoterome_receipt_sha256": sha256_file(promoterome / "receipt.json"),
            "promoter_windows_sha256": sha256_file(promoterome / "promoter_windows.tsv"),
            "promoter_transcripts_sha256": sha256_file(promoterome / "promoter_transcripts.tsv"),
            "promoter_fasta_sha256": sha256_file(promoterome / "promoter_windows.fa"),
        },
        "assembly_id": args.assembly_id, "upstream_bp": args.upstream_bp,
        "downstream_bp": args.downstream_bp,
        "sequence_orientation": "transcript_5prime_to_3prime",
        "record_policy": "One record per distinct genomic TSS window and strand; shared transcripts are listed in the FASTA header.",
        "files": manifest_files, "sha256sums_sha256": sha256_file(sums),
        "non_claims": [
            "Transcript annotation does not establish TSS usage in the assayed cells.",
            "A TSS window is sequence context, not proof of promoter activity.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                                          encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--promoterome", type=Path, required=True)
    parser.add_argument("--gene", action="append", required=True)
    parser.add_argument("--upstream-bp", type=int, default=500)
    parser.add_argument("--downstream-bp", type=int, default=200)
    parser.add_argument("--expected-genome-id", required=True)
    parser.add_argument("--assembly-id", default="GRCh38")
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = export(args)
    print(json.dumps({"status": "ok", "files": len(manifest["files"]),
                      "records": sum(item["record_count"] for item in manifest["files"]),
                      "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
