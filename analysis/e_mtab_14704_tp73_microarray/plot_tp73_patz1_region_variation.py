#!/usr/bin/env python3
"""Plot TP73/PATZ1 regional Clariom D summaries with replicate variation.

The input CSVs are produced by the local E-MTAB-14704 probe-level workflow.
They already contain per-probe-set-region group means and standard deviations
across the three CEL files per condition.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
EXPORTS = Path("/home/clawbio/.openclaw/workspace/exports")

GROUPS = [
    ("GFP", "#4b4b4b"),
    ("DNp73beta", "#2878b5"),
    ("TAp73alpha", "#c23b22"),
]
CONTRASTS = [
    ("TAp73alpha vs GFP", "log2FC_TAp73alpha_vs_GFP", "#c23b22"),
    ("DNp73beta vs GFP", "log2FC_DNp73beta_vs_GFP", "#2878b5"),
    ("TAp73alpha vs DNp73beta", "log2FC_TAp73alpha_vs_DNp73beta", "#6b4fa3"),
]


@dataclass(frozen=True)
class GenePlotSpec:
    symbol: str
    csv_name: str
    chrom: str
    strand: int
    ensembl_transcript: str
    exons: tuple[tuple[int, int], ...]


GENES = [
    GenePlotSpec(
        symbol="TP73",
        csv_name="tp73_probeset_region_intensity_chrom_order.csv",
        chrom="chr1",
        strand=1,
        ensembl_transcript="ENST00000378295",
        exons=(
            (3652516, 3652641),
            (3682333, 3682430),
            (3683060, 3683180),
            (3685525, 3685627),
            (3695469, 3695568),
            (3696904, 3697035),
            (3702417, 3702579),
            (3723256, 3723393),
            (3728589, 3728736),
            (3729327, 3729448),
            (3730000, 3730148),
            (3730927, 3731065),
            (3731463, 3731556),
            (3732747, 3736201),
        ),
    ),
    GenePlotSpec(
        symbol="PATZ1",
        csv_name="patz1_probeset_region_intensity_chrom_order.csv",
        chrom="chr22",
        strand=-1,
        ensembl_transcript="ENST00000266269",
        exons=(
            (31325804, 31327309),
            (31328787, 31328924),
            (31335692, 31335863),
            (31342897, 31342960),
            (31344332, 31346346),
        ),
    ),
]


def load_psr(spec: GenePlotSpec) -> pd.DataFrame:
    df = pd.read_csv(ROOT / spec.csv_name)
    df = df[df["feature_type"] == "probeset_region"].copy()
    df["midpoint"] = (df["start"] + df["stop"]) / 2.0
    return df.sort_values("midpoint")


def three_prime_exon_indices(spec: GenePlotSpec, n: int = 3) -> set[int]:
    ordered = sorted(enumerate(spec.exons), key=lambda item: item[1][0])
    if spec.strand == 1:
        return {idx for idx, _ in ordered[-n:]}
    return {idx for idx, _ in ordered[:n]}


def draw_exon_track(ax: plt.Axes, spec: GenePlotSpec, xmin: float, xmax: float) -> None:
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlim(xmin, xmax)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_alpha(0.25)
    ax.hlines(0.5, xmin, xmax, color="#7c7c7c", linewidth=1.0, alpha=0.7)

    highlighted = three_prime_exon_indices(spec)
    for idx, (start, stop) in enumerate(spec.exons):
        left = max(start, xmin)
        right = min(stop, xmax)
        if right < xmin or left > xmax:
            continue
        color = "#25364a" if idx in highlighted else "#b7c4d2"
        height = 0.34 if idx in highlighted else 0.24
        ax.add_patch(
            plt.Rectangle(
                (left, 0.5 - height / 2),
                max(1.0, right - left),
                height,
                facecolor=color,
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
        )

    if spec.strand == 1:
        arrow_start, arrow_end, ha = xmax - (xmax - xmin) * 0.16, xmax, "left"
    else:
        arrow_start, arrow_end, ha = xmin + (xmax - xmin) * 0.16, xmin, "right"
    ax.annotate(
        "3' end",
        xy=(arrow_end, 0.82),
        xytext=(arrow_start, 0.82),
        ha=ha,
        va="center",
        fontsize=9,
        color="#25364a",
        arrowprops={"arrowstyle": "->", "color": "#25364a", "linewidth": 1.0},
    )
    strand_label = "+" if spec.strand == 1 else "-"
    ax.text(
        xmin,
        0.08,
        f"{spec.symbol} canonical exons ({spec.ensembl_transcript}, {strand_label} strand)",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#4d5965",
    )


def draw_gene_column(fig: plt.Figure, axes: np.ndarray, spec: GenePlotSpec) -> None:
    df = load_psr(spec)
    xmin = float(min(df["start"].min(), min(start for start, _ in spec.exons)))
    xmax = float(max(df["stop"].max(), max(stop for _, stop in spec.exons)))
    span = xmax - xmin

    draw_exon_track(axes[0], spec, xmin, xmax)

    offsets = np.linspace(-0.006, 0.006, len(GROUPS)) * span
    for offset, (group, color) in zip(offsets, GROUPS, strict=True):
        axes[1].errorbar(
            df["midpoint"] + offset,
            df[f"{group}_mean_log2"],
            yerr=df[f"{group}_sd_log2"],
            fmt="none",
            ecolor=color,
            elinewidth=1.5,
            capsize=2.6,
            alpha=0.9,
            label=f"{group} mean +/- SD",
        )
        axes[1].plot(
            df["midpoint"] + offset,
            df[f"{group}_mean_log2"],
            linestyle="none",
            marker="_",
            markersize=7,
            markeredgewidth=1.7,
            color=color,
        )
    axes[1].set_ylabel("Regional log2 intensity")
    axes[1].set_title(f"{spec.symbol}: probe-set-region summaries with CEL-file variation")
    axes[1].grid(True, alpha=0.22)
    axes[1].legend(frameon=False, fontsize=8, ncol=1, loc="upper left")

    axes[2].axhline(0, color="#333333", linewidth=0.8)
    for label, column, color in CONTRASTS:
        axes[2].plot(
            df["midpoint"],
            df[column],
            linewidth=1.4,
            marker="o",
            markersize=3.2,
            label=label,
            color=color,
        )
    axes[2].set_ylabel("log2 fold change")
    axes[2].set_xlabel(f"{spec.chrom} genomic coordinate (hg38)")
    axes[2].grid(True, alpha=0.22)
    axes[2].legend(frameon=False, fontsize=8, loc="best")

    for ax in axes:
        ax.set_xlim(xmin - 0.01 * span, xmax + 0.01 * span)


def main() -> None:
    EXPORTS.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(
        3,
        2,
        figsize=(16, 8.6),
        sharex="col",
        gridspec_kw={"height_ratios": [0.55, 2.3, 1.7], "hspace": 0.16, "wspace": 0.15},
    )
    for col, spec in enumerate(GENES):
        draw_gene_column(fig, ax[:, col], spec)

    fig.suptitle(
        "Clariom D regional summaries for TP73 and PATZ1: group consensus/divergence across CEL files",
        fontsize=14,
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        "Vertical bars show SD across the three CEL files per sample group. Exon track uses Ensembl canonical transcript coordinates; dark exons mark the 3' end side.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#4d5965",
    )
    out_png = ROOT / "tp73_patz1_region_variation_3prime_exons.png"
    out_svg = ROOT / "tp73_patz1_region_variation_3prime_exons.svg"
    export_png = EXPORTS / out_png.name
    export_svg = EXPORTS / out_svg.name
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(export_png, dpi=220, bbox_inches="tight")
    fig.savefig(export_svg, bbox_inches="tight")
    print(out_png)
    print(export_png)


if __name__ == "__main__":
    main()
