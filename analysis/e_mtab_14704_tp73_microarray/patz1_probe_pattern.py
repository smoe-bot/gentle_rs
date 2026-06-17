#!/usr/bin/env python3
"""Chromosome-ordered PATZ1 probe/probeset intensity summaries."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
FEATURE_TABLE = ROOT / "all_arrays_raw_features.tsv"
SQLITE = Path(
    "/home/clawbio/.openclaw/workspace/cache/bioconductor/"
    "pd_clariom_d_human/pd.clariom.d.human/inst/extdata/"
    "pd.clariom.d.human.sqlite"
)
PATZ1_TC = "TC2200008477.hg.1"


def load_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    con = sqlite3.connect(SQLITE)
    regions = pd.read_sql_query(
        """
        select man_fsetid, fsetid, type, start, stop
        from featureSet
        where transcript_cluster_id=?
          and type in (1, 2)
        order by coalesce(start, 999999999), fsetid
        """,
        con,
        params=(PATZ1_TC,),
    )
    regions["feature_type"] = regions["type"].map({1: "junction", 2: "probeset_region"})
    probes = pd.read_sql_query(
        """
        select pm.fid as probe_id, pm.fsetid, pm.atom, pm.x, pm.y
        from pmfeature pm
        join featureSet fs on fs.fsetid = pm.fsetid
        where fs.transcript_cluster_id=?
          and fs.type in (1, 2)
        order by coalesce(fs.start, 999999999), fs.fsetid, pm.fid
        """,
        con,
        params=(PATZ1_TC,),
    )
    return regions, probes


def load_quantile_normalized_selected(selected_probe_ids: set[int]) -> tuple[list[str], pd.DataFrame]:
    df = pd.read_csv(FEATURE_TABLE, sep="\t", comment="#")
    sample_cols = [c for c in df.columns if c.endswith(".CEL")]
    probe_ids = df["probe_id"].to_numpy(dtype=np.int64)
    raw = df[sample_cols].to_numpy(dtype=np.float32, copy=True)
    del df

    sorted_sum = np.zeros(raw.shape[0], dtype=np.float64)
    for col in range(raw.shape[1]):
        sorted_sum += np.sort(raw[:, col].astype(np.float64, copy=False))
    mean_sorted = (sorted_sum / raw.shape[1]).astype(np.float32)
    del sorted_sum

    selected_sorted = np.array(sorted(selected_probe_ids), dtype=np.int64)
    probe_order = np.argsort(probe_ids)
    sorted_probe_ids = probe_ids[probe_order]
    positions = np.searchsorted(sorted_probe_ids, selected_sorted)
    ok = (positions < len(sorted_probe_ids)) & (sorted_probe_ids[positions] == selected_sorted)
    selected_sorted = selected_sorted[ok]
    raw_positions = probe_order[positions[ok]]

    values = np.empty((len(selected_sorted), raw.shape[1]), dtype=np.float32)
    for col in range(raw.shape[1]):
        order = np.argsort(raw[:, col], kind="mergesort")
        ranks = np.empty(raw.shape[0], dtype=np.int32)
        ranks[order] = np.arange(raw.shape[0], dtype=np.int32)
        values[:, col] = mean_sorted[ranks[raw_positions]]
        del ranks, order

    out = pd.DataFrame(values, columns=sample_cols)
    out.insert(0, "probe_id", selected_sorted)
    return sample_cols, out


def add_contrasts(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    groups = {
        "DNp73beta": [c for c in sample_cols if "AdDNp73beta" in c],
        "GFP": [c for c in sample_cols if "AdGFP" in c],
        "TAp73alpha": [c for c in sample_cols if "AdTAp73alpha" in c],
    }
    out = df.copy()
    for col in sample_cols:
        out[f"log2_{col}"] = np.log2(out[col].astype(float) + 1.0)
    for group, cols in groups.items():
        log_cols = [f"log2_{c}" for c in cols]
        out[f"{group}_mean_log2"] = out[log_cols].mean(axis=1)
        out[f"{group}_sd_log2"] = out[log_cols].std(axis=1)
    out["log2FC_TAp73alpha_vs_GFP"] = out["TAp73alpha_mean_log2"] - out["GFP_mean_log2"]
    out["log2FC_DNp73beta_vs_GFP"] = out["DNp73beta_mean_log2"] - out["GFP_mean_log2"]
    out["log2FC_TAp73alpha_vs_DNp73beta"] = (
        out["TAp73alpha_mean_log2"] - out["DNp73beta_mean_log2"]
    )
    return out


def main() -> None:
    regions, probes = load_regions()
    sample_cols, probe_values = load_quantile_normalized_selected(set(probes["probe_id"]))

    probe_table = probes.merge(regions, on="fsetid", how="left").merge(
        probe_values, on="probe_id", how="left"
    )
    probe_table = add_contrasts(probe_table, sample_cols)
    probe_table = probe_table.sort_values(["start", "stop", "man_fsetid", "probe_id"])
    probe_table.to_csv(ROOT / "patz1_pm_probe_intensity_chrom_order.csv", index=False)

    psr_rows = []
    for _, region in regions.iterrows():
        sub = probe_table[probe_table["fsetid"] == region["fsetid"]]
        if sub.empty:
            continue
        row = region.to_dict()
        row["n_pm_probes"] = len(sub)
        for col in sample_cols:
            row[col] = float(np.nanmedian(sub[col]))
        psr_rows.append(row)
    psr = pd.DataFrame(psr_rows)
    psr = add_contrasts(psr, sample_cols).sort_values(["start", "stop", "man_fsetid"])
    psr.to_csv(ROOT / "patz1_probeset_region_intensity_chrom_order.csv", index=False)

    plot_df = psr[psr["feature_type"] == "probeset_region"].copy()
    plot_df["midpoint"] = (plot_df["start"] + plot_df["stop"]) / 2
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for group, color in [("GFP", "#555555"), ("DNp73beta", "#2878b5"), ("TAp73alpha", "#c23b22")]:
        ax1.plot(
            plot_df["midpoint"],
            plot_df[f"{group}_mean_log2"],
            marker="o",
            linewidth=1.5,
            label=group,
            color=color,
        )
    ax1.set_ylabel("Mean log2 intensity")
    ax1.set_title("PATZ1 Clariom D probe-set regions, chr22 order")
    ax1.legend(frameon=False)
    ax1.grid(True, alpha=0.25)

    ax2.axhline(0, color="#333333", linewidth=0.8)
    ax2.plot(
        plot_df["midpoint"],
        plot_df["log2FC_TAp73alpha_vs_GFP"],
        marker="o",
        label="TAp73alpha vs GFP",
        color="#c23b22",
    )
    ax2.plot(
        plot_df["midpoint"],
        plot_df["log2FC_DNp73beta_vs_GFP"],
        marker="o",
        label="DNp73beta vs GFP",
        color="#2878b5",
    )
    ax2.plot(
        plot_df["midpoint"],
        plot_df["log2FC_TAp73alpha_vs_DNp73beta"],
        marker="o",
        label="TAp73alpha vs DNp73beta",
        color="#6b4fa3",
    )
    ax2.set_xlabel("chr22 genomic coordinate (hg38)")
    ax2.set_ylabel("log2 fold change")
    ax2.legend(frameon=False)
    ax2.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "patz1_probe_region_intensity_chrom_order.png", dpi=180)
    fig.savefig(ROOT / "patz1_probe_region_intensity_chrom_order.svg")

    print(f"Wrote {len(probe_table)} PM-probe rows and {len(psr)} feature-region rows.")


if __name__ == "__main__":
    main()
