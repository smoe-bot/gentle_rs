#!/usr/bin/env python3
"""Prepare five consistent 30-track locus-report requests from validated templates."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys

try:
    from .prepare_regulatory_region_indexes import require
except ImportError:
    from prepare_regulatory_region_indexes import require


TARGET_TRANSCRIPTS = {
    "CD44": ["ENST00000428726", "ENST00000263398", "ENST00000278386"],
    "TGFB1": ["ENST00000221930", "ENST00001090430", "ENST00001114525"],
    "SERPINE1": ["ENST00000223095", "ENST00000870828", "ENST00000950058"],
    "PATZ1": ["ENST00000930048", "ENST00000405309", "ENST00000266269"],
    "TP73": ["ENST00000378295", "ENST00000378288", "ENST00000378290"],
}

MATRICES = [
    ("TP53", "MA0106.3"), ("TP63", "MA0525.2"), ("TP73", "MA0861.2"),
    ("SP1", "MA0079.5"), ("CGGBP1", "MA2538.1"), ("HELT", "MA2628.1"),
    ("IRF9", "MA0653.1"), ("PROX1", "MA0794.1"), ("FOXF2", "MA0030.2"),
    ("MYBL2", "MA0777.1"), ("PLAG1", "MA0163.1"), ("STAT1", "MA0137.4"),
    ("PATZ1", "MA1961.2"), ("ELK4", "MA0076.3"), ("REST", "MA0138.3"),
    ("LMX1B", "MA0703.3"), ("POU2F2", "MA0507.3"),
    ("TFAP2C", "MA0524.3"), ("TFAP2C", "MA0814.3"), ("TFAP2C", "MA0815.1"),
    ("ZBED1", "MA0749.2"), ("TGIF2", "MA0797.1"), ("GATA5", "MA0766.3"),
    ("HES7", "MA0822.1"), ("KLF15", "MA1513.2"), ("NFKB2", "MA0778.2"),
    ("POU2F1", "MA0785.2"), ("SATB1", "MA1963.2"), ("IRF2", "MA0051.2"),
    ("MEF2B", "MA0660.1"),
]

COLORS = ["#7b2cbf", "#ca6702", "#005f73", "#ae2012", "#1b4332", "#2a9d8f",
          "#6a4c93", "#8a5a44", "#e63946", "#3f37c9", "#457b9d", "#2b9348",
          "#bc6c25", "#0077b6", "#d00000", "#606c38", "#dda15e", "#9c6644",
          "#b56576", "#6d597a", "#118ab2", "#588157", "#e76f51", "#8338ec",
          "#3a86ff", "#ff006e", "#4361ee", "#7209b7", "#0081a7", "#f07167"]


def score_tracks() -> list[dict]:
    tracks = []
    for index, ((factor, matrix), color) in enumerate(zip(MATRICES, COLORS, strict=True), 1):
        suffix = matrix.lower().replace(".", "_")
        label = f"{factor} JASPAR PWM"
        if factor == "TFAP2C":
            label += f" ({matrix})"
        tracks.append({
            "track_id": f"{factor.lower()}_{suffix}_jaspar_2026", "label": label,
            "provider_kind": "jaspar_pwm", "source_ids": [matrix],
            "factors": [{"factor_id": factor, "factor_label": factor}],
            "score_kind": "llr_background_tail_log10", "calibration_state": "matrix_specific",
            "calibration_statement": "JASPAR 2026 CORE matrix-specific PWM score; not calibrated as biochemical affinity.",
            "strand_policy": "both", "clip_negative": True, "display_threshold": 0,
            "top_hit_count": 5, "scale_mode": "independent", "color_hint": color,
            "display_order": index,
        })
    return tracks


def prepare(args: argparse.Namespace) -> None:
    template_dir = args.template_dir.resolve(strict=True)
    entries = args.entry_dir.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "output directory must be absent or empty")
    output.mkdir(parents=True, exist_ok=True)
    outputs = output / "outputs"
    outputs.mkdir()
    base = json.loads((template_dir / "CD44.json").read_text())
    tracks = score_tracks()
    require(len(tracks) == 30 and len({t["track_id"] for t in tracks}) == 30,
            "30-track panel must have unique track IDs")
    for gene, transcripts in TARGET_TRANSCRIPTS.items():
        template_path = template_dir / f"{gene}.json"
        request = (json.loads(template_path.read_text()) if template_path.exists()
                   else copy.deepcopy(base))
        entry_path = (entries / f"{gene}.ensembl116.entry.json").resolve(strict=True)
        entry = json.loads(entry_path.read_text())
        present = {row["transcript_id"] for row in entry["transcripts"]}
        require(set(transcripts) == present, f"offline entry transcript selection mismatch: {gene}")
        stem = f"{gene}_luciferase_planning_EnsemblReg_30TF"
        request.update({
            "gene_query": gene, "allow_ensembl_network": False,
            "ensembl_entry_path": str(entry_path),
            "annotation_release": "Ensembl 116 prepared GRCh38; selected reporter transcripts",
            "output_panel_id": f"{gene}_panel", "output_seq_id": f"{gene}_locus",
            "svg_path": str((outputs / f"{stem}.svg").resolve()),
            "png_path": str((outputs / f"{stem}.png").resolve()),
            "pdf_path": str((outputs / f"{stem}.pdf").resolve()),
            "receipt_path": str((outputs / f"{stem}.receipt.json").resolve()),
            "display_report_path": str((outputs / f"{stem}.report.json").resolve()),
            "reporter_report_path": str((outputs / f"{stem}.reporter.json").resolve()),
            "regulatory_score_tracks": tracks,
        })
        reporter = copy.deepcopy(request["reporter_architecture_request"])
        reporter.update({"seq_id": "", "gene_label": gene, "transcript_ids": transcripts})
        request["reporter_architecture_request"] = reporter
        (output / f"{gene}.json").write_text(json.dumps(request, indent=2) + "\n")
    (output / "target_transcripts.json").write_text(
        json.dumps(TARGET_TRANSCRIPTS, indent=2, sort_keys=True) + "\n")
    (output / "jaspar_30_track_panel.json").write_text(
        json.dumps({"schema": "gentle.jaspar_target_panel.v1", "tracks": tracks}, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-dir", type=Path, required=True)
    parser.add_argument("--entry-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args)
    print(json.dumps({"status": "ok", "targets": len(TARGET_TRANSCRIPTS),
                      "tracks": len(MATRICES), "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
