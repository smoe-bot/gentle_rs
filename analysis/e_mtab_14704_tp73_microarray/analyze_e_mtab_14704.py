#!/usr/bin/env python3
"""First-pass targeted analysis of E-MTAB-14704 Clariom D arrays.

This script uses the local Bioconductor pd.clariom.d.human platform package
and the APT feature-intensity dump produced from the nine CEL files.
"""

from __future__ import annotations

import math
import re
import sqlite3
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
FEATURE_TABLE = ROOT / "all_arrays_raw_features.tsv"
SQLITE = Path(
    "/home/clawbio/.openclaw/workspace/cache/bioconductor/"
    "pd_clariom_d_human/pd.clariom.d.human/inst/extdata/"
    "pd.clariom.d.human.sqlite"
)
NETAFFX_TRANSCRIPT = SQLITE.with_name("netaffxTranscript.rda")

LIST_A = [
    "TFAP2E",
    "TFAP2B",
    "KLF15",
    "TFAP2A",
    "TFAP2C",
    "PLAGL2",
    "EGR3",
    "NFIX",
    "HES7",
    "PATZ1",
    "ELK3",
    "GLIS3",
    "INSM1",
    "E2F1",
    "SP1",
]
LIST_B = ["CDX4", "ISL1", "POU4F2", "POU2F2", "BCL6B", "FOXF2", "POU1F1", "LMX1B"]
FUS_PROJECT = ["FUS", "HDAC1", "HDAC2", "HDAC6", "TERT", "REST", "MDM2", "TARDBP", "E2F1"]
ALL_TARGETS = sorted(set(LIST_A + LIST_B + FUS_PROJECT + ["TP73"]))

TYPE_MAP = {1: "junction", 2: "probeset_region", 3: "unmapped"}


@dataclass(frozen=True)
class FeatureSet:
    symbol: str
    list_name: str
    transcript_cluster_id: str
    feature_id: int
    feature_label: str
    feature_type: str
    seqname: str
    strand: str
    start: float
    stop: float
    n_probes: int


def split_gene_symbols(geneassignment: object) -> set[str]:
    symbols: set[str] = set()
    if geneassignment is None or (isinstance(geneassignment, float) and math.isnan(geneassignment)):
        return symbols
    for record in str(geneassignment).split(" /// "):
        fields = [part.strip() for part in record.split(" // ")]
        if len(fields) >= 2 and fields[1] and fields[1] != "---":
            symbols.add(fields[1])
    return symbols


def target_list_name(symbol: str) -> str:
    names = []
    if symbol in LIST_A:
        names.append("List A")
    if symbol in LIST_B:
        names.append("List B")
    if symbol in FUS_PROJECT:
        names.append("TP73-FUS panel")
    if symbol == "TP73":
        names.append("TP73")
    return "; ".join(names)


def load_transcript_annotation() -> pd.DataFrame:
    cache = ROOT / "netaffx_transcript_targets.csv"
    if cache.exists():
        return pd.read_csv(cache)

    warnings.filterwarnings("ignore", category=UserWarning, module="rdata")
    obj = rdata.conversion.convert(rdata.parser.parse_file(NETAFFX_TRANSCRIPT))["netaffxTranscript"]
    df = obj.data.copy()
    rows = []
    for _, row in df.iterrows():
        symbols = split_gene_symbols(row["geneassignment"])
        for symbol in sorted(symbols & set(ALL_TARGETS)):
            rows.append(
                {
                    "symbol": symbol,
                    "list_name": target_list_name(symbol),
                    "transcript_cluster_id": row["transcriptclusterid"],
                    "probesetid": row["probesetid"],
                    "seqname": row["seqname"],
                    "strand": row["strand"],
                    "start": row["start"],
                    "stop": row["stop"],
                    "totalprobes": row["totalprobes"],
                    "category": row["category"],
                    "geneassignment": row["geneassignment"],
                }
            )
    out = pd.DataFrame(rows).sort_values(["symbol", "transcript_cluster_id"])
    out.to_csv(cache, index=False)
    return out


def fetch_probe_ids(con: sqlite3.Connection, fsetid: int) -> list[int]:
    rows = con.execute("select fid from pmfeature where fsetid=? order by fid", (int(fsetid),)).fetchall()
    return [int(r[0]) for r in rows]


def core_features(con: sqlite3.Connection, targets: pd.DataFrame) -> tuple[list[FeatureSet], dict[int, list[int]]]:
    features: list[FeatureSet] = []
    probes_by_feature: dict[int, list[int]] = {}
    for row in targets.itertuples(index=False):
        core_rows = con.execute(
            "select fsetid from core_mps where transcript_cluster_id=?",
            (row.transcript_cluster_id,),
        ).fetchall()
        for (fsetid,) in core_rows:
            probe_ids = fetch_probe_ids(con, int(fsetid))
            if not probe_ids:
                continue
            probes_by_feature[int(fsetid)] = probe_ids
            features.append(
                FeatureSet(
                    symbol=row.symbol,
                    list_name=row.list_name,
                    transcript_cluster_id=row.transcript_cluster_id,
                    feature_id=int(fsetid),
                    feature_label=row.transcript_cluster_id,
                    feature_type="transcript_cluster_core",
                    seqname=row.seqname,
                    strand=row.strand,
                    start=float(row.start),
                    stop=float(row.stop),
                    n_probes=len(probe_ids),
                )
            )
    return features, probes_by_feature


def tp73_probe_region_features(con: sqlite3.Connection) -> tuple[list[FeatureSet], dict[int, list[int]]]:
    rows = con.execute(
        """
        select man_fsetid, fsetid, type, start, stop
        from featureSet
        where transcript_cluster_id='TC0100006620.hg.1'
          and type in (1, 2)
        order by coalesce(start, 999999999), fsetid
        """
    ).fetchall()
    features: list[FeatureSet] = []
    probes_by_feature: dict[int, list[int]] = {}
    for label, fsetid, type_id, start, stop in rows:
        probe_ids = fetch_probe_ids(con, int(fsetid))
        if not probe_ids:
            continue
        probes_by_feature[int(fsetid)] = probe_ids
        features.append(
            FeatureSet(
                symbol="TP73",
                list_name="TP73 probe pattern",
                transcript_cluster_id="TC0100006620.hg.1",
                feature_id=int(fsetid),
                feature_label=label,
                feature_type=TYPE_MAP.get(type_id, f"type_{type_id}"),
                seqname="chr1",
                strand="+",
                start=float(start) if start is not None else np.nan,
                stop=float(stop) if stop is not None else np.nan,
                n_probes=len(probe_ids),
            )
        )
    return features, probes_by_feature


def load_quantile_normalized_selected(selected_probe_ids: set[int]) -> tuple[list[str], dict[int, np.ndarray]]:
    print("Reading raw feature matrix...")
    df = pd.read_csv(FEATURE_TABLE, sep="\t", comment="#")
    sample_cols = [c for c in df.columns if c.endswith(".CEL")]
    probe_ids = df["probe_id"].to_numpy(dtype=np.int64)
    raw = df[sample_cols].to_numpy(dtype=np.float32, copy=True)
    del df

    print("Quantile-normalizing selected probes...")
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
    selected_position_set = dict(zip(raw_positions.tolist(), selected_sorted.tolist(), strict=True))

    norm_by_probe: dict[int, np.ndarray] = {
        int(pid): np.empty(raw.shape[1], dtype=np.float32) for pid in selected_sorted
    }
    for col in range(raw.shape[1]):
        order = np.argsort(raw[:, col], kind="mergesort")
        ranks = np.empty(raw.shape[0], dtype=np.int32)
        ranks[order] = np.arange(raw.shape[0], dtype=np.int32)
        for pos, pid in selected_position_set.items():
            norm_by_probe[int(pid)][col] = mean_sorted[ranks[pos]]
        del ranks, order

    return sample_cols, norm_by_probe


def welch(values: np.ndarray, left: list[int], right: list[int]) -> float:
    if len(left) < 2 or len(right) < 2:
        return np.nan
    p = stats.ttest_ind(values[left], values[right], equal_var=False, nan_policy="omit").pvalue
    return float(p) if np.isfinite(p) else np.nan


def paired_t(values: np.ndarray, left: list[int], right: list[int]) -> float:
    """Paired t-test using replicate suffix order, if both groups have n>=2."""
    if len(left) != len(right) or len(left) < 2:
        return np.nan
    p = stats.ttest_rel(values[left], values[right], nan_policy="omit").pvalue
    return float(p) if np.isfinite(p) else np.nan


def summarize_features(
    features: list[FeatureSet],
    probes_by_feature: dict[int, list[int]],
    norm_by_probe: dict[int, np.ndarray],
    sample_cols: list[str],
) -> pd.DataFrame:
    groups = {
        "DNp73beta": [i for i, c in enumerate(sample_cols) if "AdDNp73beta" in c],
        "GFP": [i for i, c in enumerate(sample_cols) if "AdGFP" in c],
        "TAp73alpha": [i for i, c in enumerate(sample_cols) if "AdTAp73alpha" in c],
    }
    rows = []
    for feature in features:
        matrix = [norm_by_probe[pid] for pid in probes_by_feature[feature.feature_id] if pid in norm_by_probe]
        if not matrix:
            continue
        expr = np.log2(np.vstack(matrix).astype(np.float64) + 1.0)
        values = np.nanmedian(expr, axis=0)
        row = {
            "symbol": feature.symbol,
            "list_name": feature.list_name,
            "transcript_cluster_id": feature.transcript_cluster_id,
            "feature_label": feature.feature_label,
            "feature_type": feature.feature_type,
            "seqname": feature.seqname,
            "strand": feature.strand,
            "start": feature.start,
            "stop": feature.stop,
            "n_probes": feature.n_probes,
        }
        for group, idx in groups.items():
            row[f"{group}_mean_log2"] = float(np.nanmean(values[idx]))
            row[f"{group}_sd_log2"] = float(np.nanstd(values[idx], ddof=1))
        row["log2FC_TAp73alpha_vs_GFP"] = row["TAp73alpha_mean_log2"] - row["GFP_mean_log2"]
        row["log2FC_DNp73beta_vs_GFP"] = row["DNp73beta_mean_log2"] - row["GFP_mean_log2"]
        row["log2FC_TAp73alpha_vs_DNp73beta"] = row["TAp73alpha_mean_log2"] - row["DNp73beta_mean_log2"]
        row["p_TAp73alpha_vs_GFP"] = welch(values, groups["TAp73alpha"], groups["GFP"])
        row["p_DNp73beta_vs_GFP"] = welch(values, groups["DNp73beta"], groups["GFP"])
        row["p_TAp73alpha_vs_DNp73beta"] = welch(values, groups["TAp73alpha"], groups["DNp73beta"])
        row["p_paired_TAp73alpha_vs_GFP"] = paired_t(values, groups["TAp73alpha"], groups["GFP"])
        row["p_paired_DNp73beta_vs_GFP"] = paired_t(values, groups["DNp73beta"], groups["GFP"])
        row["p_paired_TAp73alpha_vs_DNp73beta"] = paired_t(
            values, groups["TAp73alpha"], groups["DNp73beta"]
        )
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in [
        "p_TAp73alpha_vs_GFP",
        "p_DNp73beta_vs_GFP",
        "p_TAp73alpha_vs_DNp73beta",
        "p_paired_TAp73alpha_vs_GFP",
        "p_paired_DNp73beta_vs_GFP",
        "p_paired_TAp73alpha_vs_DNp73beta",
    ]:
        if col in out and out[col].notna().any():
            mask = out[col].notna()
            out.loc[mask, "q_" + col[2:]] = multipletests(out.loc[mask, col], method="fdr_bh")[1]
    return out


def gene_rollup(core_df: pd.DataFrame) -> pd.DataFrame:
    numeric = [
        "DNp73beta_mean_log2",
        "GFP_mean_log2",
        "TAp73alpha_mean_log2",
        "log2FC_TAp73alpha_vs_GFP",
        "log2FC_DNp73beta_vs_GFP",
        "log2FC_TAp73alpha_vs_DNp73beta",
        "p_TAp73alpha_vs_GFP",
        "p_DNp73beta_vs_GFP",
        "p_TAp73alpha_vs_DNp73beta",
        "p_paired_TAp73alpha_vs_GFP",
        "p_paired_DNp73beta_vs_GFP",
        "p_paired_TAp73alpha_vs_DNp73beta",
    ]
    rows = []
    for symbol, sub in core_df.groupby("symbol", sort=True):
        preferred = sub.copy()
        # Prefer ordinary coding/main transcript clusters over circRNA companion probes.
        if symbol == "E2F1":
            coding = preferred[preferred["transcript_cluster_id"] == "TC2000008894.hg.1"]
            if not coding.empty:
                preferred = coding
        row = {
            "symbol": symbol,
            "list_name": "; ".join(sorted(set("; ".join(preferred["list_name"]).split("; ")))),
            "n_transcript_clusters": len(preferred),
            "transcript_cluster_ids": ";".join(preferred["transcript_cluster_id"]),
        }
        for col in numeric:
            row[col] = float(np.nanmedian(preferred[col]))
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in [
        "p_TAp73alpha_vs_GFP",
        "p_DNp73beta_vs_GFP",
        "p_TAp73alpha_vs_DNp73beta",
        "p_paired_TAp73alpha_vs_GFP",
        "p_paired_DNp73beta_vs_GFP",
        "p_paired_TAp73alpha_vs_DNp73beta",
    ]:
        mask = out[col].notna()
        out.loc[mask, "q_" + col[2:]] = multipletests(out.loc[mask, col], method="fdr_bh")[1]
    return out.sort_values(["list_name", "symbol"])


def main() -> None:
    con = sqlite3.connect(SQLITE)
    targets = load_transcript_annotation()
    targets.to_csv(ROOT / "target_transcript_clusters.csv", index=False)

    core, core_probes = core_features(con, targets)
    tp73_pattern, tp73_probes = tp73_probe_region_features(con)
    all_probes = set()
    for probe_list in list(core_probes.values()) + list(tp73_probes.values()):
        all_probes.update(probe_list)
    print(f"Selected {len(all_probes)} unique probes from {len(core)} core and {len(tp73_pattern)} TP73 regional features.")

    sample_cols, norm_by_probe = load_quantile_normalized_selected(all_probes)
    pd.Series(sample_cols, name="sample").to_csv(ROOT / "samples.csv", index=False)

    core_df = summarize_features(core, core_probes, norm_by_probe, sample_cols)
    tp73_df = summarize_features(tp73_pattern, tp73_probes, norm_by_probe, sample_cols)

    gene_df = gene_rollup(core_df)
    list_ab_df = gene_df[gene_df["symbol"].isin(set(LIST_A + LIST_B))].copy()
    fus_df = gene_df[gene_df["symbol"].isin(FUS_PROJECT)].copy()

    core_df.to_csv(ROOT / "target_transcript_cluster_summary.csv", index=False)
    gene_df.to_csv(ROOT / "target_gene_summary.csv", index=False)
    list_ab_df.to_csv(ROOT / "list_a_b_gene_summary.csv", index=False)
    fus_df.to_csv(ROOT / "tp73_fus_panel_gene_summary.csv", index=False)
    tp73_df.to_csv(ROOT / "tp73_regional_probe_pattern.csv", index=False)

    psr = tp73_df[tp73_df["feature_type"] == "probeset_region"].copy()
    psr["max_abs_variant_log2FC_vs_GFP"] = psr[
        ["log2FC_TAp73alpha_vs_GFP", "log2FC_DNp73beta_vs_GFP"]
    ].abs().max(axis=1)
    psr.sort_values("max_abs_variant_log2FC_vs_GFP", ascending=False).head(20).to_csv(
        ROOT / "tp73_top_changing_probe_regions.csv", index=False
    )

    print("Wrote:")
    for path in [
        "target_transcript_clusters.csv",
        "target_transcript_cluster_summary.csv",
        "target_gene_summary.csv",
        "list_a_b_gene_summary.csv",
        "tp73_fus_panel_gene_summary.csv",
        "tp73_regional_probe_pattern.csv",
        "tp73_top_changing_probe_regions.csv",
    ]:
        print(f"  {ROOT / path}")


if __name__ == "__main__":
    main()
