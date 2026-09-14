#!/usr/bin/env python3
"""Build a checksum-bound Ensembl/secondary-annotation TSS comparison.

The primary annotation is read from an existing TSS-profile report.  The
secondary annotation is read from a GENtle prepared-genome transcript index.
No transcript equivalence or consensus TSS is inferred.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import html
import json
from pathlib import Path
from typing import Any


REPORT_SCHEMA = "gentle.transcript_start_annotation_comparison.v1"
RECEIPT_SCHEMA = "gentle.transcript_start_annotation_comparison_receipt.v1"
TSS_REPORT_SCHEMA = "gentle.tss_tfbs_profiles.v1"
PAGE_WIDTH = 1400
PLOT_LEFT = 255
PLOT_RIGHT = 1050


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_bytes())
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def chromosome_token(raw: str) -> str:
    value = raw.strip()
    lower = value.lower().removeprefix("chr")
    accession = lower.split(".", 1)[0]
    if accession.startswith("nc_") and accession[3:].isdigit():
        number = int(accession[3:])
        return {23: "x", 24: "y", 12920: "mt"}.get(number, str(number))
    if lower in {"m", "mt"}:
        return "mt"
    return lower.lstrip("0") or "0"


def transcript_tss(row: dict[str, Any]) -> int:
    strand = row.get("strand")
    require(strand in {"+", "-"}, "secondary transcript lacks a +/- strand")
    key = "transcript_start_1based" if strand == "+" else "transcript_end_1based"
    value = row.get(key)
    require(isinstance(value, int) and value > 0, "secondary transcript has invalid geometry")
    return value


def grouped_secondary(
    rows: list[dict[str, Any]], gene: str, chromosome: str, strand: str
) -> list[dict[str, Any]]:
    selected = []
    for row in rows:
        if str(row.get("gene_name", "")).casefold() != gene.casefold():
            continue
        if chromosome_token(str(row.get("chromosome", ""))) != chromosome_token(chromosome):
            continue
        if row.get("strand") != strand:
            continue
        selected.append(row)
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in selected:
        transcript_id = row.get("transcript_id")
        require(isinstance(transcript_id, str) and transcript_id, "secondary transcript ID missing")
        grouped[transcript_tss(row)].append(transcript_id)
    return [
        {"tss_1based": tss, "transcript_ids": sorted(set(ids)), "transcript_count": len(set(ids))}
        for tss, ids in sorted(grouped.items())
    ]


def signed_transcript_delta(primary_tss: int, secondary_tss: int, strand: str) -> int:
    # Positive means the secondary start lies downstream in transcript direction.
    return secondary_tss - primary_tss if strand == "+" else primary_tss - secondary_tss


def nearest_cluster(tss: int, clusters: list[dict[str, Any]], strand: str) -> dict[str, Any] | None:
    if not clusters:
        return None
    nearest = min(clusters, key=lambda row: (abs(row["tss_1based"] - tss), row["tss_1based"]))
    delta = signed_transcript_delta(tss, nearest["tss_1based"], strand)
    return {
        "tss_1based": nearest["tss_1based"],
        "delta_bp_transcript_direction": delta,
        "absolute_delta_bp": abs(delta),
        "relationship": "exact" if delta == 0 else "nearest",
    }


def catalog_identity(catalog: dict[str, Any], genome_id: str) -> dict[str, Any]:
    row = catalog.get(genome_id)
    require(isinstance(row, dict), f"catalog has no object for {genome_id}")
    identity = {"genome_id": genome_id, "description": row.get("description")}
    for key in (
        "ncbi_taxonomy_id", "ncbi_assembly_accession", "ncbi_assembly_name",
        "sequence_remote", "annotations_remote", "ensembl_template",
    ):
        if key in row:
            identity[key] = row[key]
    return identity


def assembly_core(identity: dict[str, Any]) -> str | None:
    name = identity.get("ncbi_assembly_name")
    if isinstance(name, str) and name:
        return name.split(".", 1)[0].casefold()
    template = identity.get("ensembl_template")
    if isinstance(template, dict):
        stem = template.get("file_stem")
        if isinstance(stem, str) and "." in stem:
            return stem.rsplit(".", 1)[-1].casefold()
    return None


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    primary = read_object(args.primary_tss_report)
    require(primary.get("schema") == TSS_REPORT_SCHEMA, "unexpected primary TSS report schema")
    catalog = read_object(args.catalog)
    manifest = read_object(args.secondary_manifest)
    secondary_id = manifest.get("genome_id")
    require(isinstance(secondary_id, str) and secondary_id, "secondary manifest genome_id missing")
    primary_id = primary.get("reference", {}).get("genome_id")
    require(isinstance(primary_id, str) and primary_id, "primary report genome_id missing")
    primary_identity = catalog_identity(catalog, primary_id)
    secondary_identity = catalog_identity(catalog, secondary_id)
    primary_assembly = str(primary.get("reference", {}).get("assembly", "")).casefold()
    secondary_assembly = assembly_core(secondary_identity)
    require(primary_assembly and secondary_assembly == primary_assembly,
            "primary and secondary annotations are not on the same assembly core")
    transcript_path = Path(str(manifest.get("transcript_index_path", "")))
    if not transcript_path.is_absolute():
        transcript_path = args.secondary_manifest.parent / transcript_path
    require(transcript_path.is_file(), "secondary transcript index is missing")
    annotation_path = Path(str(manifest.get("annotation_path", "")))
    if not annotation_path.is_absolute():
        annotation_path = args.secondary_manifest.parent / annotation_path
    require(annotation_path.is_file(), "secondary annotation is missing")
    declared_annotation_sha1 = manifest.get("annotation_sha1")
    require(isinstance(declared_annotation_sha1, str)
            and declared_annotation_sha1 == sha1(annotation_path),
            "secondary annotation SHA-1 does not match its GENtle manifest")
    transcripts = json.loads(transcript_path.read_bytes())
    require(isinstance(transcripts, list), "secondary transcript index is not an array")

    requested = list(dict.fromkeys(args.gene))
    require(requested, "at least one --gene is required")
    gene_reports = []
    windows = primary.get("windows")
    require(isinstance(windows, list), "primary report windows must be an array")
    for gene in requested:
        gene_windows = [w for w in windows if w.get("record", {}).get("gene_symbol") == gene]
        require(gene_windows, f"primary report contains no windows for {gene}")
        strands = {w["record"]["geometry"]["strand"] for w in gene_windows}
        chromosomes = {w["record"]["geometry"]["chromosome"] for w in gene_windows}
        require(len(strands) == 1 and len(chromosomes) == 1,
                f"primary {gene} windows do not share chromosome/strand")
        strand = strands.pop()
        chromosome = chromosomes.pop()
        primary_clusters = []
        primary_transcript_ids: set[str] = set()
        for window in sorted(gene_windows, key=lambda w: w["record"]["geometry"]["tss_1based"]):
            record = window["record"]
            ids = record.get("transcripts", [])
            require(isinstance(ids, list) and ids, f"primary {gene} window lacks transcripts")
            primary_transcript_ids.update(ids)
            primary_clusters.append({
                "promoter_id": record["promoter_id"],
                "tss_1based": record["geometry"]["tss_1based"],
                "transcript_ids": ids,
                "transcript_count": len(ids),
                "selected_for_reporter": bool(window.get("selected")),
            })
        secondary_clusters = grouped_secondary(transcripts, gene, chromosome, strand)
        require(secondary_clusters, f"secondary annotation contains no {gene} transcript clusters")
        for cluster in primary_clusters:
            cluster["nearest_secondary"] = nearest_cluster(
                cluster["tss_1based"], secondary_clusters, strand
            )
        for cluster in secondary_clusters:
            cluster["nearest_primary"] = nearest_cluster(
                cluster["tss_1based"], primary_clusters, strand
            )
        gene_reports.append({
            "gene_symbol": gene,
            "chromosome": chromosome,
            "strand": strand,
            "primary_transcript_count": len(primary_transcript_ids),
            "secondary_transcript_count": sum(row["transcript_count"] for row in secondary_clusters),
            "primary_tss_cluster_count": len(primary_clusters),
            "secondary_tss_cluster_count": len(secondary_clusters),
            "exact_shared_tss_count": len(
                {row["tss_1based"] for row in primary_clusters}
                & {row["tss_1based"] for row in secondary_clusters}
            ),
            "primary_clusters": primary_clusters,
            "secondary_clusters": secondary_clusters,
        })

    return {
        "schema": REPORT_SCHEMA,
        "coordinate_system": "1-based closed genomic coordinates",
        "comparison_rule": (
            "TSS is transcript start: genomic start on + and genomic end on -. "
            "Nearest relationships are coordinate descriptions, not transcript equivalence."
        ),
        "primary": primary_identity,
        "secondary": secondary_identity,
        "inputs": {
            "primary_tss_report": {"path": str(args.primary_tss_report), "sha256": sha256(args.primary_tss_report)},
            "catalog": {"path": str(args.catalog), "sha256": sha256(args.catalog)},
            "secondary_manifest": {"path": str(args.secondary_manifest), "sha256": sha256(args.secondary_manifest)},
            "secondary_annotation": {
                "path": str(annotation_path),
                "sha1": declared_annotation_sha1,
                "sha256": sha256(annotation_path),
                "source": manifest.get("annotation_source"),
            },
            "secondary_transcript_index": {"path": str(transcript_path), "sha256": sha256(transcript_path)},
        },
        "genes": gene_reports,
        "producer": {
            "revision": args.producer_revision,
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "gentle_cli": str(args.gentle_cli.resolve()),
            "gentle_cli_sha256": sha256(args.gentle_cli.resolve()),
        },
        "non_claims": (
            "Annotation agreement or proximity does not establish TSS usage, transcript abundance, "
            "promoter activity, CUT&RUN occupancy, TF binding, or reporter function."
        ),
    }


def marker_x(tss: int, start: int, end: int) -> float:
    if end == start:
        return (PLOT_LEFT + PLOT_RIGHT) / 2
    return PLOT_LEFT + (tss - start) * (PLOT_RIGHT - PLOT_LEFT) / (end - start)


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def render_gene_svg(report: dict[str, Any], gene: dict[str, Any]) -> str:
    all_tss = [r["tss_1based"] for r in gene["primary_clusters"] + gene["secondary_clusters"]]
    pad = max(100, (max(all_tss) - min(all_tss)) // 20)
    axis_start = min(all_tss) - pad
    axis_end = max(all_tss) + pad
    selected = [r for r in gene["primary_clusters"] if r["selected_for_reporter"]]
    table_rows = max(len(selected), len(gene["secondary_clusters"]))
    height = 500 + 24 * table_rows
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{PAGE_WIDTH}" height="{height}" '
        f'viewBox="0 0 {PAGE_WIDTH} {height}" data-gentle-schema="{REPORT_SCHEMA}" '
        f'data-gentle-gene="{esc(gene["gene_symbol"])}" data-gentle-plot-left="{PLOT_LEFT}" '
        f'data-gentle-plot-right="{PLOT_RIGHT}">',
        '<style>text{font-family:"DejaVu Sans",sans-serif;fill:#243238}.title{font-size:24px;font-weight:700}'
        '.sub{font-size:12px;fill:#526268}.head{font-size:14px;font-weight:700;fill:#173f43}'
        '.lab{font-size:11px}.small{font-size:9px;fill:#596d72}.axis{stroke:#60767a;stroke-width:1}'
        '.ens{stroke:#2f6ea3;stroke-width:2}.sel{stroke:#db6b3f;stroke-width:4}.ref{stroke:#4d8b55;stroke-width:2}'
        '.grid{stroke:#d9e2df;stroke-width:1}</style>',
        '<rect width="100%" height="100%" fill="#fbfaf4"/>',
        '<rect x="0" y="0" width="100%" height="94" fill="#eaf1ed"/>',
        f'<text class="title" x="36" y="38">{esc(gene["gene_symbol"])} transcript-start annotation comparison</text>',
        f'<text class="sub" x="36" y="64">chr{esc(gene["chromosome"])} · strand {esc(gene["strand"])} · shared genomic axis · no consensus TSS inferred</text>',
        f'<text class="sub" x="36" y="82">Primary: {esc(report["primary"]["genome_id"])} | Secondary: {esc(report["secondary"]["genome_id"])}</text>',
        f'<text class="head" x="36" y="125">Inventory</text>',
        f'<text class="lab" x="36" y="149">Ensembl: {gene["primary_transcript_count"]} transcripts, {gene["primary_tss_cluster_count"]} TSS coordinates</text>',
        f'<text class="lab" x="440" y="149">RefSeq: {gene["secondary_transcript_count"]} transcripts, {gene["secondary_tss_cluster_count"]} TSS coordinates</text>',
        f'<text class="lab" x="840" y="149">Exact shared coordinates: {gene["exact_shared_tss_count"]}</text>',
        f'<line class="axis" x1="{PLOT_LEFT}" x2="{PLOT_RIGHT}" y1="214" y2="214"/>',
        f'<text class="small" x="{PLOT_LEFT}" y="201">{axis_start:,}</text>',
        f'<text class="small" x="{PLOT_RIGHT}" y="201" text-anchor="end">{axis_end:,}</text>',
        '<text class="lab" x="36" y="236">Ensembl</text>',
        '<text class="lab" x="36" y="282">NCBI RefSeq</text>',
    ]
    for row in gene["primary_clusters"]:
        x = marker_x(row["tss_1based"], axis_start, axis_end)
        cls = "sel" if row["selected_for_reporter"] else "ens"
        out.append(f'<line class="{cls}" x1="{x:.3f}" x2="{x:.3f}" y1="220" y2="248" data-tss="{row["tss_1based"]}"/>')
    for row in gene["secondary_clusters"]:
        x = marker_x(row["tss_1based"], axis_start, axis_end)
        out.append(f'<line class="ref" x1="{x:.3f}" x2="{x:.3f}" y1="266" y2="294" data-tss="{row["tss_1based"]}"/>')
    out.extend([
        '<text class="small" x="1080" y="236">orange = reporter-selected Ensembl TSS</text>',
        '<text class="head" x="36" y="334">Reporter-selected Ensembl starts and nearest RefSeq start</text>',
        '<text class="small" x="36" y="354">Ensembl TSS</text><text class="small" x="190" y="354">transcripts</text>'
        '<text class="small" x="370" y="354">nearest RefSeq TSS</text><text class="small" x="560" y="354">transcript-direction delta</text>',
    ])
    y = 378
    for row in selected:
        nearest = row["nearest_secondary"]
        delta = nearest["delta_bp_transcript_direction"]
        delta_label = "exact" if delta == 0 else f'{delta:+d} bp'
        out.extend([
            f'<line class="grid" x1="36" x2="1040" y1="{y + 6}" y2="{y + 6}"/>',
            f'<text class="lab" x="36" y="{y}">{row["tss_1based"]:,}</text>',
            f'<text class="lab" x="190" y="{y}">{row["transcript_count"]}</text>',
            f'<text class="lab" x="370" y="{y}">{nearest["tss_1based"]:,}</text>',
            f'<text class="lab" x="560" y="{y}">{delta_label}</text>',
            f'<text class="small" x="690" y="{y}">{esc(row["promoter_id"])}</text>',
        ])
        y += 24
    secondary_x = 1080
    out.append(f'<text class="head" x="{secondary_x}" y="334">RefSeq TSS inventory</text>')
    for idx, row in enumerate(gene["secondary_clusters"]):
        nearest = row["nearest_primary"]
        delta = nearest["delta_bp_transcript_direction"]
        label = "exact" if delta == 0 else f'{delta:+d} bp to nearest Ensembl'
        out.append(f'<text class="small" x="{secondary_x}" y="{378 + idx * 24}">{row["tss_1based"]:,} · n={row["transcript_count"]} · {label}</text>')
    footer_y = height - 52
    out.extend([
        f'<text class="small" x="36" y="{footer_y}">Positive delta means the RefSeq start lies downstream in transcript direction; negative means upstream.</text>',
        f'<text class="small" x="36" y="{footer_y + 18}">{esc(report["non_claims"])}</text>',
        '</svg>',
    ])
    return "".join(out)


def write_outputs(args: argparse.Namespace, report: dict[str, Any]) -> None:
    require(not args.output_report.exists(), f"output already exists: {args.output_report}")
    require(not args.output_receipt.exists(), f"output already exists: {args.output_receipt}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.output_report.parent.mkdir(parents=True, exist_ok=True)
    args.output_report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs: dict[str, str] = {args.output_report.name: sha256(args.output_report)}
    pages = []
    for gene in report["genes"]:
        path = args.output_dir / f'{gene["gene_symbol"]}_Ensembl_RefSeq_TSS_comparison.svg'
        require(not path.exists(), f"output already exists: {path}")
        path.write_text(render_gene_svg(report, gene), encoding="utf-8")
        outputs[path.name] = sha256(path)
        pages.append({"gene_symbol": gene["gene_symbol"], "path": str(path), "sha256": outputs[path.name]})
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "report_schema": REPORT_SCHEMA,
        "report_sha256": outputs[args.output_report.name],
        "pages": pages,
        "outputs": outputs,
        "producer": report["producer"],
    }
    args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-tss-report", type=Path, required=True)
    parser.add_argument("--secondary-manifest", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--gene", action="append", required=True)
    parser.add_argument("--producer-revision", required=True)
    parser.add_argument("--gentle-cli", type=Path, required=True)
    parser.add_argument("--output-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-receipt", type=Path, required=True)
    args = parser.parse_args()
    write_outputs(args, build_report(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
