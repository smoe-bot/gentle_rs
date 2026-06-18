#!/usr/bin/env python3
"""Build a compact Clariom D probe-set activity presentation.

The E-MTAB-14704 CEL-derived intermediates are intentionally not committed as
bulk data. This script derives an inspectable, reproducible view from explicit
local inputs: an APT raw PM-probe feature table, the Clariom D SQLite platform
annotation, and the vendor probeset annotation ZIP.
"""

from __future__ import annotations

import csv
import html
import io
import json
import math
import sqlite3
import zipfile
import argparse
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "analysis" / "e_mtab_14704_tp73_microarray"
ANNOTATION = (
    ROOT
    / "test_files"
    / "fixtures"
    / "affymetrix_clariom_d_human_na36_hg38_subset"
    / "clariom_d_human_na36_hg38_gene_panel.probesets.tsv"
)
VENDOR_PROBESET_ZIP = (
    ROOT
    / "data"
    / "publication_resources"
    / "rostock_p73_clariomd_e_mtab_14704"
    / "library"
    / "Clariom_D_Human-na36-hg38-probeset-csv.zip"
)
VENDOR_PROBESET_MEMBER = "Clariom_D_Human.na36.hg38.probeset.csv"
SQLITE = WORK / "Rlib" / "pd.clariom.d.human" / "extdata" / "pd.clariom.d.human.sqlite"
RAW_FEATURES = WORK / "all_arrays_raw_features.tsv"
OUT = WORK / "gene_panel_probe_set_activity"

GENES = ["TP73", "FUS", "PATZ1", "E2F1", "TARDBP", "PLK1", "TERT", "HDAC1", "HDAC2", "HDAC6"]
GROUPS = {
    "DNp73beta": [
        "P_SKMel29_AdDNp73beta_1.CEL",
        "P_SKMel29_AdDNp73beta_2.CEL",
        "P_SKMel29_AdDNp73beta_3.CEL",
    ],
    "GFP": [
        "P_SKMel29_AdGFP_1.CEL",
        "P_SKMel29_AdGFP_2.CEL",
        "P_SKMel29_AdGFP_3.CEL",
    ],
    "TAp73alpha": [
        "P_SKMel29_AdTAp73alpha_1.CEL",
        "P_SKMel29_AdTAp73alpha_2.CEL",
        "P_SKMel29_AdTAp73alpha_3.CEL",
    ],
}
SAMPLES = [sample for samples in GROUPS.values() for sample in samples]
CONTRAST_COLUMNS = [
    "log2_TAp73alpha_minus_GFP",
    "log2_DNp73beta_minus_GFP",
    "log2_TAp73alpha_minus_DNp73beta",
]
PAIR_REPLICATES = ("1", "2", "3")
PROBESET_COLUMNS = [
    "gene_symbols",
    "probeset_id",
    "seqname",
    "strand",
    "start",
    "stop",
    "probe_count",
    "transcript_cluster_id",
    "locus_type",
    "exon_id",
    "psr_id",
    "probeset_type",
    "psr_type",
    "junction_start_edge",
    "junction_stop_edge",
    "level",
    "has_cds",
]


def clean(value: str | None) -> str:
    if value is None or value == "---":
        return ""
    return value


def split_fixture_symbols(value: str | None) -> set[str]:
    symbols: set[str] = set()
    for chunk in (value or "").replace(" /// ", ";").split(";"):
        chunk = chunk.strip()
        if chunk and chunk != "---":
            symbols.add(chunk)
    return symbols


def split_vendor_assignment(value: str | None) -> set[str]:
    symbols: set[str] = set()
    for item in (value or "").split(" /// "):
        parts = [part.strip() for part in item.split(" // ")]
        if len(parts) > 1 and parts[1] and parts[1] != "---":
            symbols.add(parts[1])
    return symbols


def parse_float(value: str) -> float:
    return float(value) if value else math.nan


def mean(values: list[float]) -> float:
    observed = [v for v in values if not math.isnan(v)]
    return sum(observed) / len(observed) if observed else math.nan


def log2(value: float) -> float:
    return math.log2(value + 1.0) if not math.isnan(value) else math.nan


def fmt(value: float) -> str:
    return "" if math.isnan(value) else f"{value:.4f}"


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def read_probe_sets() -> list[dict[str, str]]:
    wanted = set(GENES)
    rows: list[dict[str, str]] = []
    if ANNOTATION.exists():
        with ANNOTATION.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                symbols = split_fixture_symbols(row["gene_symbols"]) & wanted
                if not symbols:
                    continue
                out = {column: row.get(column, "") for column in PROBESET_COLUMNS}
                out["gene_symbols"] = ";".join(sorted(symbols, key=GENES.index))
                rows.append(out)
    observed = set().union(*(split_fixture_symbols(row["gene_symbols"]) for row in rows)) if rows else set()
    missing = wanted - observed
    if missing:
        with zipfile.ZipFile(VENDOR_PROBESET_ZIP) as archive:
            with archive.open(VENDOR_PROBESET_MEMBER) as zipped:
                text = io.TextIOWrapper(zipped, encoding="utf-8-sig", errors="replace", newline="")
                reader = csv.DictReader(line for line in text if not line.startswith("#"))
                for row in reader:
                    symbols = split_vendor_assignment(row.get("gene_assignment")) & missing
                    if not symbols:
                        continue
                    out = {column: row.get(column, "") for column in PROBESET_COLUMNS}
                    out["gene_symbols"] = ";".join(sorted(symbols, key=GENES.index))
                    rows.append(out)
    order = {gene: i for i, gene in enumerate(GENES)}
    return sorted(
        rows,
        key=lambda row: (
            min(order[symbol] for symbol in split_fixture_symbols(row["gene_symbols"])),
            clean(row.get("start")) or "999999999999",
            row["probeset_id"],
        ),
    )


def fetch_pm_features(probe_sets: list[dict[str, str]]) -> dict[str, list[dict[str, int]]]:
    wanted = [row["probeset_id"] for row in probe_sets]
    placeholders = ",".join("?" for _ in wanted)
    query = f"""
        SELECT fs.man_fsetid, pm.fid, pm.x, pm.y
        FROM featureSet fs
        JOIN pmfeature pm ON pm.fsetid = fs.fsetid
        WHERE fs.man_fsetid IN ({placeholders})
        ORDER BY fs.man_fsetid, pm.fid
    """
    conn = sqlite3.connect(SQLITE)
    try:
        by_probe_set: dict[str, list[dict[str, int]]] = defaultdict(list)
        for probeset_id, fid, x, y in conn.execute(query, wanted):
            by_probe_set[probeset_id].append({"probe_id": int(fid), "x": int(x), "y": int(y)})
        return by_probe_set
    finally:
        conn.close()


def read_raw_intensities(probe_ids: set[int]) -> dict[int, dict[str, float]]:
    intensities: dict[int, dict[str, float]] = {}
    with RAW_FEATURES.open(newline="") as handle:
        header: list[str] | None = None
        for line in handle:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            raise RuntimeError(f"No header found in {RAW_FEATURES}")
        index = {name: i for i, name in enumerate(header)}
        missing_samples = [sample for sample in SAMPLES if sample not in index]
        if missing_samples:
            raise RuntimeError(f"Missing sample columns: {missing_samples}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            probe_id = int(fields[index["probe_id"]])
            if probe_id not in probe_ids:
                continue
            intensities[probe_id] = {sample: parse_float(fields[index[sample]]) for sample in SAMPLES}
            if len(intensities) == len(probe_ids):
                break
    return intensities


def write_tsv(path: Path, rows: list[dict[str, str | int | float]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_probe_sets(
    probe_sets: list[dict[str, str]],
    pm_features: dict[str, list[dict[str, int]]],
    intensities: dict[int, dict[str, float]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    probe_level_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    for row in probe_sets:
        probeset_id = row["probeset_id"]
        gene_symbols = split_fixture_symbols(row["gene_symbols"])
        primary_gene = min(gene_symbols, key=GENES.index)
        features = pm_features.get(probeset_id, [])
        observed = [feature for feature in features if feature["probe_id"] in intensities]
        for feature in observed:
            values = intensities[feature["probe_id"]]
            probe_level_rows.append(
                {
                    "gene": row["gene_symbols"],
                    "primary_gene": primary_gene,
                    "probeset_id": probeset_id,
                    "probe_id": str(feature["probe_id"]),
                    "x": str(feature["x"]),
                    "y": str(feature["y"]),
                    **{sample: fmt(values[sample]) for sample in SAMPLES},
                }
            )
        sample_means = {
            sample: mean([intensities[feature["probe_id"]][sample] for feature in observed])
            for sample in SAMPLES
        }
        group_log2 = {
            group: mean([log2(sample_means[sample]) for sample in samples])
            for group, samples in GROUPS.items()
        }
        summary_rows.append(
            {
                "gene": row["gene_symbols"],
                "primary_gene": primary_gene,
                "probeset_id": probeset_id,
                "probeset_type": row["probeset_type"],
                "seqname": row["seqname"],
                "strand": row["strand"],
                "start": clean(row["start"]),
                "stop": clean(row["stop"]),
                "transcript_cluster_id": row["transcript_cluster_id"],
                "exon_id": clean(row["exon_id"]),
                "junction_start_edge": clean(row["junction_start_edge"]),
                "junction_stop_edge": clean(row["junction_stop_edge"]),
                "annotated_probe_count": row["probe_count"],
                "pm_probe_count": str(len(features)),
                "observed_pm_probe_count": str(len(observed)),
                **{sample: fmt(sample_means[sample]) for sample in SAMPLES},
                "log2_mean_DNp73beta": fmt(group_log2["DNp73beta"]),
                "log2_mean_GFP": fmt(group_log2["GFP"]),
                "log2_mean_TAp73alpha": fmt(group_log2["TAp73alpha"]),
                "log2_TAp73alpha_minus_GFP": fmt(group_log2["TAp73alpha"] - group_log2["GFP"]),
                "log2_DNp73beta_minus_GFP": fmt(group_log2["DNp73beta"] - group_log2["GFP"]),
                "log2_TAp73alpha_minus_DNp73beta": fmt(
                    group_log2["TAp73alpha"] - group_log2["DNp73beta"]
                ),
            }
        )
    return summary_rows, probe_level_rows


def colour(value: float, lo: float, hi: float) -> str:
    if math.isnan(value):
        return "#f3f4f6"
    if hi <= lo:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - lo) / (hi - lo)))
    # Blue-white-red without a single-hue palette.
    if t < 0.5:
        f = t / 0.5
        r = int(49 + (245 - 49) * f)
        g = int(104 + (245 - 104) * f)
        b = int(176 + (245 - 176) * f)
    else:
        f = (t - 0.5) / 0.5
        r = int(245 + (190 - 245) * f)
        g = int(245 + (55 - 245) * f)
        b = int(245 + (45 - 245) * f)
    return f"#{r:02x}{g:02x}{b:02x}"


def write_svg(rows: list[dict[str, str]]) -> None:
    row_h = 14
    label_w = 292
    cell_w = 26
    top = 88
    columns = SAMPLES + CONTRAST_COLUMNS
    values = [
        parse_float(row[column])
        for row in rows
        for column in columns
        if row.get(column)
    ]
    lo, hi = min(values), max(values)
    width = label_w + cell_w * len(columns) + 220
    height = top + row_h * len(rows) + 84
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">',
        "<title>TP73 overexpression probe-set raw PM activity</title>",
        "<style>"
        "text{font-family:Inter,Arial,sans-serif;font-size:10px;fill:#111827}"
        ".small{font-size:9px;fill:#374151}.tiny{font-size:8px;fill:#4b5563}"
        ".gene{font-weight:700}.axis{stroke:#d1d5db;stroke-width:1}"
        "</style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="18" y="24" style="font-size:16px;font-weight:700">Probe-set activity in TP73 isoform overexpression arrays</text>',
        '<text x="18" y="43" class="small">Raw PM-probe intensities summarized by Clariom D probeset; geometry/triage view, not an isoform-support call.</text>',
    ]
    for i, column in enumerate(columns):
        x = label_w + i * cell_w + 18
        parts.append(
            f'<text x="{x}" y="78" class="tiny" transform="rotate(-45 {x} 78)">'
            f"{html.escape(column.replace('P_SKMel29_Ad', '').replace('.CEL', ''))}</text>"
        )
    last_gene = None
    for r, row in enumerate(rows):
        y = top + r * row_h
        gene = row["gene"]
        if gene != last_gene:
            parts.append(f'<line x1="16" y1="{y - 5}" x2="{width - 18}" y2="{y - 5}" class="axis"/>')
            parts.append(f'<text x="18" y="{y + 9}" class="gene">{html.escape(gene)}</text>')
            last_gene = gene
        label = f"{row['probeset_id']} {row['probeset_type']} {row['start']}-{row['stop']}".strip()
        parts.append(f'<text x="66" y="{y + 9}" class="tiny">{html.escape(label[:56])}</text>')
        for c, column in enumerate(columns):
            x = label_w + c * cell_w
            value = parse_float(row[column]) if row[column] else math.nan
            fill = colour(value, lo, hi)
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w - 1}" height="{row_h - 1}" '
                f'fill="{fill}"><title>{html.escape(row["probeset_id"])} {column}: {fmt(value)}</title></rect>'
            )
    legend_x = label_w + cell_w * len(columns) + 34
    parts.extend(
        [
            f'<text x="{legend_x}" y="92" class="small">Scale spans raw means and log2 contrasts</text>',
            f'<text x="{legend_x}" y="110" class="tiny">min {lo:.2f}</text>',
            f'<text x="{legend_x + 92}" y="110" class="tiny">max {hi:.2f}</text>',
        ]
    )
    for i in range(80):
        v = lo + (hi - lo) * i / 79
        parts.append(
            f'<rect x="{legend_x + i}" y="116" width="1" height="12" fill="{colour(v, lo, hi)}"/>'
        )
    parts.append("</svg>")
    (OUT / "probe_set_activity_heatmap.svg").write_text("\n".join(parts) + "\n")


def write_html(rows: list[dict[str, str]], meta: dict[str, object]) -> None:
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_gene[row["gene"]].append(row)
    cards = []
    for gene in GENES:
        gene_rows = [row for row in rows if row["primary_gene"] == gene]
        top_ta = sorted(
            gene_rows,
            key=lambda row: parse_float(row["log2_TAp73alpha_minus_GFP"]),
            reverse=True,
        )[:5]
        cards.append(f"<section><h2>{html.escape(gene)}</h2>")
        cards.append(
            f"<p>{len(gene_rows)} probesets, "
            f"{sum(int(r['observed_pm_probe_count']) for r in gene_rows)} observed PM probes.</p>"
        )
        cards.append("<table><thead><tr><th>probeset</th><th>type</th><th>TA-GFP</th><th>DN-GFP</th><th>locus</th></tr></thead><tbody>")
        for row in top_ta:
            locus = f"{row['seqname']}:{row['start']}-{row['stop']}" if row["start"] else "junction"
            cards.append(
                "<tr>"
                f"<td>{html.escape(row['probeset_id'])}</td>"
                f"<td>{html.escape(row['probeset_type'])}</td>"
                f"<td>{html.escape(row['log2_TAp73alpha_minus_GFP'])}</td>"
                f"<td>{html.escape(row['log2_DNp73beta_minus_GFP'])}</td>"
                f"<td>{html.escape(locus)}</td>"
                "</tr>"
            )
        cards.append("</tbody></table></section>")
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TP73 overexpression probe-set activity</title>
<style>
body {{ font-family: Inter, Arial, sans-serif; margin: 28px; color: #111827; background: #fff; }}
h1 {{ font-size: 24px; margin-bottom: 6px; }}
h2 {{ font-size: 17px; margin: 22px 0 6px; }}
p {{ max-width: 980px; color: #374151; }}
table {{ border-collapse: collapse; width: 100%; max-width: 1120px; font-size: 13px; }}
th, td {{ border-bottom: 1px solid #e5e7eb; padding: 5px 7px; text-align: left; }}
th {{ background: #f9fafb; }}
.figure {{ margin: 18px 0; border: 1px solid #e5e7eb; overflow-x: auto; }}
.meta {{ font-size: 12px; color: #4b5563; }}
</style>
</head>
<body>
<h1>Probe-set activity in TP73 isoform overexpression arrays</h1>
<p>This is a figure-preparation surface for raw PM-probe activity summarized by Clariom D probeset across three GFP controls, three DNp73beta arrays, and three TAp73alpha arrays. It is a conservative geometry/activity view, not an isoform-support call.</p>
<p class="meta">Generated from {html.escape(display_path(RAW_FEATURES))} and {html.escape(display_path(SQLITE))}. {html.escape(json.dumps(meta, sort_keys=True))}</p>
<div class="figure"><img src="probe_set_activity_heatmap.svg" alt="Probe-set activity heatmap"></div>
{''.join(cards)}
</body>
</html>
"""
    (OUT / "index.html").write_text(html_text)


def paired_log2(row: dict[str, str], replicate: str, condition: str) -> float:
    control = log2(parse_float(row[f"P_SKMel29_AdGFP_{replicate}.CEL"]))
    treated = log2(parse_float(row[f"P_SKMel29_Ad{condition}_{replicate}.CEL"]))
    return treated - control


def write_paired_gene_summary(rows: list[dict[str, str]]) -> None:
    fieldnames = ["primary_gene", "probesets", "pm_probes"]
    for replicate in PAIR_REPLICATES:
        fieldnames.extend(
            [
                f"E{replicate}_median_DN_minus_GFP",
                f"E{replicate}_median_TA_minus_GFP",
            ]
        )
    fieldnames.extend(
        [
            "median_DN_minus_GFP_across_pairs",
            "median_TA_minus_GFP_across_pairs",
        ]
    )
    with (OUT / "paired_gene_level_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in GENES:
            gene_rows = [row for row in rows if row["primary_gene"] == gene]
            if not gene_rows:
                continue
            out: dict[str, str | int] = {
                "primary_gene": gene,
                "probesets": len(gene_rows),
                "pm_probes": sum(int(row["observed_pm_probe_count"]) for row in gene_rows),
            }
            dn_medians: list[float] = []
            ta_medians: list[float] = []
            for replicate in PAIR_REPLICATES:
                dn = sorted(paired_log2(row, replicate, "DNp73beta") for row in gene_rows)
                ta = sorted(paired_log2(row, replicate, "TAp73alpha") for row in gene_rows)
                dn_median = dn[len(dn) // 2]
                ta_median = ta[len(ta) // 2]
                out[f"E{replicate}_median_DN_minus_GFP"] = fmt(dn_median)
                out[f"E{replicate}_median_TA_minus_GFP"] = fmt(ta_median)
                dn_medians.append(dn_median)
                ta_medians.append(ta_median)
            out["median_DN_minus_GFP_across_pairs"] = fmt(sorted(dn_medians)[len(dn_medians) // 2])
            out["median_TA_minus_GFP_across_pairs"] = fmt(sorted(ta_medians)[len(ta_medians) // 2])
            writer.writerow(out)


def write_matplotlib_figures(rows: list[dict[str, str]]) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.colors import TwoSlopeNorm
    except Exception as exc:  # pragma: no cover - optional plotting dependency.
        print(f"skipping PNG/PDF plots; matplotlib unavailable: {exc}")
        return []

    written: list[str] = []
    ordered: list[dict[str, str]] = []
    for gene in GENES:
        gene_rows = [row for row in rows if row["primary_gene"] == gene]
        gene_rows.sort(
            key=lambda row: sum(paired_log2(row, replicate, "TAp73alpha") for replicate in PAIR_REPLICATES)
            / len(PAIR_REPLICATES),
            reverse=True,
        )
        ordered.extend(gene_rows)

    # Gene-level contrast distribution.
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    for ax, (column, title, colour_name) in zip(
        axes,
        [
            ("log2_TAp73alpha_minus_GFP", "TAp73alpha vs GFP", "#b91c1c"),
            ("log2_DNp73beta_minus_GFP", "DNp73beta vs GFP", "#1d4ed8"),
        ],
    ):
        for idx, gene in enumerate(GENES):
            values = [parse_float(row[column]) for row in rows if row["primary_gene"] == gene]
            xs = [idx + ((n % 11) - 5) * 0.025 for n in range(len(values))]
            ax.scatter(xs, values, s=18, c=colour_name, alpha=0.55, edgecolors="none")
            if values:
                median = sorted(values)[len(values) // 2]
                ax.plot([idx - 0.3, idx + 0.3], [median, median], color="black", linewidth=2)
                ax.text(
                    idx,
                    max(values) + (0.35 if max(values) > 2 else 0.22),
                    f"n={len(values)}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="#374151",
                )
        ax.axhline(0, color="#6b7280", linewidth=1)
        ax.axhline(2, color="#dc2626", linewidth=1, linestyle="--")
        ax.axhline(-2, color="#2563eb", linewidth=1, linestyle="--")
        ax.set_ylabel("log2 contrast")
        ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
        ax.set_ylim(-2.5, 8.1)
    axes[-1].set_xticks(range(len(GENES)))
    axes[-1].set_xticklabels(GENES, rotation=35, ha="right")
    fig.suptitle("Probe-set activity by gene: raw PM-probe means across 9 arrays", fontsize=15, fontweight="bold")
    fig.text(
        0.01,
        0.01,
        "Each point is one Clariom D probeset; black segment is median. Dashed lines mark +/-2 log2. Visual inspection only.",
        fontsize=9,
        color="#374151",
    )
    for extension in ("png", "svg", "pdf"):
        output = OUT / f"gene_contrast_probe_set_summary.{extension}"
        fig.savefig(output, dpi=180 if extension == "png" else None)
        written.append(output.name)
    plt.close(fig)

    height = max(10, len(ordered) * 0.045)
    sample_columns: list[tuple[str, str]] = []
    for replicate in PAIR_REPLICATES:
        sample_columns.extend(
            [
                (f"P_SKMel29_AdGFP_{replicate}.CEL", f"E{replicate} GFP"),
                (f"P_SKMel29_AdDNp73beta_{replicate}.CEL", f"E{replicate} DN"),
                (f"P_SKMel29_AdTAp73alpha_{replicate}.CEL", f"E{replicate} TA"),
            ]
        )
    matrix = np.array([[log2(parse_float(row[column])) for column, _ in sample_columns] for row in ordered])
    fig, ax = plt.subplots(figsize=(9.2, height), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(sample_columns)))
    ax.set_xticklabels([label for _, label in sample_columns], rotation=45, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title("Individual array probe-set activity: 10-gene panel", loc="left", fontsize=14, fontweight="bold")
    add_gene_separators(ax, ordered, text_x=-0.8, line_color="white")
    for x in (2.5, 5.5):
        ax.axvline(x, color="white", linewidth=1.2)
    colourbar = fig.colorbar(image, ax=ax, shrink=0.8, pad=0.03)
    colourbar.set_label("log2 raw PM-probe mean")
    fig.text(
        0.02,
        0.005,
        "Columns are grouped by paired experiment/time/person: GFP, DNp73beta, TAp73alpha for E1-E3.",
        fontsize=8,
    )
    for extension in ("png", "pdf"):
        output = OUT / f"probe_set_individual_arrays_heatmap_10_gene.{extension}"
        fig.savefig(output, dpi=220 if extension == "png" else None)
        written.append(output.name)
    plt.close(fig)

    paired_matrix = np.array(
        [
            [
                value
                for replicate in PAIR_REPLICATES
                for value in (
                    paired_log2(row, replicate, "DNp73beta"),
                    paired_log2(row, replicate, "TAp73alpha"),
                )
            ]
            for row in ordered
        ]
    )
    paired_labels = [
        label
        for replicate in PAIR_REPLICATES
        for label in (f"E{replicate} DN-GFP", f"E{replicate} TA-GFP")
    ]
    fig, ax = plt.subplots(figsize=(8.5, height), constrained_layout=True)
    image = ax.imshow(
        paired_matrix,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=7.5),
    )
    ax.set_xticks(range(len(paired_labels)))
    ax.set_xticklabels(paired_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title("Paired within-experiment probe-set contrasts", loc="left", fontsize=14, fontweight="bold")
    add_gene_separators(ax, ordered, text_x=-0.7, line_color="black")
    for x in (1.5, 3.5):
        ax.axvline(x, color="black", linewidth=0.8)
    colourbar = fig.colorbar(image, ax=ax, shrink=0.8, pad=0.03)
    colourbar.set_label("paired log2 contrast")
    fig.text(
        0.02,
        0.005,
        "Each contrast is computed within replicate/experiment: log2(condition PM mean + 1) - log2(matched GFP PM mean + 1).",
        fontsize=8,
    )
    for extension in ("png", "pdf"):
        output = OUT / f"probe_set_paired_contrast_heatmap_10_gene.{extension}"
        fig.savefig(output, dpi=220 if extension == "png" else None)
        written.append(output.name)
    plt.close(fig)
    return written


def add_gene_separators(ax, ordered_rows: list[dict[str, str]], text_x: float, line_color: str) -> None:
    start = 0
    for gene in GENES:
        count = sum(1 for row in ordered_rows if row["primary_gene"] == gene)
        if not count:
            continue
        ax.axhline(start - 0.5, color=line_color, linewidth=0.7)
        ax.text(text_x, start + count / 2 - 0.5, f"{gene} ({count})", va="center", ha="right", fontsize=8)
        start += count
    ax.axhline(start - 0.5, color=line_color, linewidth=0.7)


def write_readme(meta: dict[str, object]) -> None:
    text = f"""# TP73 Overexpression Probe-Set Activity

Generated local presentation for {", ".join(GENES)}.

Inputs:

- `{display_path(RAW_FEATURES)}`: APT raw PM-probe feature intensities for 9 CEL files.
- `{display_path(SQLITE)}`: local pd.clariom.d.human SQLite annotation used only to map probesets to PM probe feature ids.
- `{display_path(ANNOTATION)}`: committed small gene-panel probeset subset for genes already present there.
- `{display_path(VENDOR_PROBESET_ZIP)}`: vendor probeset annotation ZIP used to resolve genes missing from the compact fixture.

Outputs:

- `index.html`: compact presentation page.
- `probe_set_activity_heatmap.svg`: figure-ready heatmap overview.
- `probe_set_activity_summary.tsv`: probeset-level mean raw intensity and log2 group contrasts.
- `probe_level_activity.tsv`: selected PM-probe raw intensities.
- `gene_contrast_probe_set_summary.png/.svg/.pdf`: compact per-gene contrast distribution.
- `probe_set_individual_arrays_heatmap_10_gene.png/.pdf`: individual-array heatmap ordered by paired experiment.
- `probe_set_paired_contrast_heatmap_10_gene.png/.pdf`: within-experiment paired contrast heatmap.
- `paired_gene_level_summary.tsv`: per-gene medians for the paired contrasts.
- `manifest.json`: machine-readable provenance.

Caveat: the activity values are raw PM-probe intensities summarized by probeset. They are useful for visual inspection and figure preparation, but are not a formal expression model and should not be read as an isoform-support claim.

Metadata:

```json
{json.dumps(meta, indent=2, sort_keys=True)}
```
"""
    (OUT / "README.md").write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-features",
        type=Path,
        default=RAW_FEATURES,
        help="APT raw PM-probe feature table with probe_id/x/y and one column per CEL file.",
    )
    parser.add_argument(
        "--sqlite",
        type=Path,
        default=SQLITE,
        help="pd.clariom.d.human SQLite database containing featureSet and pmfeature.",
    )
    parser.add_argument(
        "--fixture-probesets",
        type=Path,
        default=ANNOTATION,
        help="Optional compact GENtle probeset fixture to use before falling back to the vendor ZIP.",
    )
    parser.add_argument(
        "--vendor-probeset-zip",
        type=Path,
        default=VENDOR_PROBESET_ZIP,
        help="Clariom D Human na36 hg38 probeset CSV ZIP for genes missing from the compact fixture.",
    )
    parser.add_argument("--output-dir", type=Path, default=OUT, help="Directory for generated local outputs.")
    parser.add_argument(
        "--genes",
        default=",".join(GENES),
        help="Comma-separated gene symbols to render, in display order.",
    )
    return parser.parse_args()


def configure_from_args(args: argparse.Namespace) -> None:
    global RAW_FEATURES, SQLITE, ANNOTATION, VENDOR_PROBESET_ZIP, OUT, GENES
    RAW_FEATURES = args.raw_features
    SQLITE = args.sqlite
    ANNOTATION = args.fixture_probesets
    VENDOR_PROBESET_ZIP = args.vendor_probeset_zip
    OUT = args.output_dir
    GENES = [gene.strip() for gene in args.genes.split(",") if gene.strip()]


def main() -> None:
    configure_from_args(parse_args())
    OUT.mkdir(parents=True, exist_ok=True)
    probe_sets = read_probe_sets()
    pm_features = fetch_pm_features(probe_sets)
    probe_ids = {feature["probe_id"] for features in pm_features.values() for feature in features}
    intensities = read_raw_intensities(probe_ids)
    missing = sorted(probe_ids - set(intensities))
    if missing:
        raise RuntimeError(f"Missing {len(missing)} PM probe ids in raw feature table")
    summary_rows, probe_level_rows = summarize_probe_sets(probe_sets, pm_features, intensities)
    summary_fields = [
        "gene",
        "primary_gene",
        "probeset_id",
        "probeset_type",
        "seqname",
        "strand",
        "start",
        "stop",
        "transcript_cluster_id",
        "exon_id",
        "junction_start_edge",
        "junction_stop_edge",
        "annotated_probe_count",
        "pm_probe_count",
        "observed_pm_probe_count",
        *SAMPLES,
        "log2_mean_DNp73beta",
        "log2_mean_GFP",
        "log2_mean_TAp73alpha",
        *CONTRAST_COLUMNS,
    ]
    probe_fields = ["gene", "primary_gene", "probeset_id", "probe_id", "x", "y", *SAMPLES]
    write_tsv(OUT / "probe_set_activity_summary.tsv", summary_rows, summary_fields)
    write_tsv(OUT / "probe_level_activity.tsv", probe_level_rows, probe_fields)
    write_svg(summary_rows)
    write_paired_gene_summary(summary_rows)
    figure_files = write_matplotlib_figures(summary_rows)
    meta = {
        "genes": GENES,
        "samples": SAMPLES,
        "probe_sets": len(summary_rows),
        "pm_probes": len(probe_level_rows),
        "groups": GROUPS,
        "paired_replicates": PAIR_REPLICATES,
        "figure_files": figure_files,
    }
    (OUT / "manifest.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    write_html(summary_rows, meta)
    write_readme(meta)
    print(f"wrote {OUT}")
    print(json.dumps(meta, sort_keys=True))


if __name__ == "__main__":
    main()
