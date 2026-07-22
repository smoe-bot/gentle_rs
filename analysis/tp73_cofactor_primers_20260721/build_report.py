#!/usr/bin/env python3
"""Summarize GENtle transcript-panel and cDNA PCR outputs for five genes."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENES = ("PATZ1", "TFAP2C", "AHCTF1", "POU2F2", "LMX1B")


def read_json(path: Path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def exon_label(hit: dict, segments: list[dict]) -> str:
    start = hit["start_0based"]
    end = hit["end_0based_exclusive"]
    ordinals = [
        str(segment["exon_ordinal"])
        for segment in segments
        if segment["local_end_0based_exclusive"] > start
        and segment["local_start_0based"] < end
    ]
    return "+".join(ordinals) if ordinals else "?"


def source_ranges(hit: dict) -> str:
    return ",".join(
        f"{row['start_0based'] + 1}-{row['end_0based_exclusive']}"
        for row in hit.get("source_ranges_0based", [])
    )


primer_rows: list[dict] = []
binding_rows: list[dict] = []
summary_rows: list[dict] = []

for gene in GENES:
    panel = read_json(ROOT / f"{gene}_endpoint_matrix.json")
    summary_rows.append(
        {
            "gene": gene,
            "completion_status": panel["completion_status"],
            "annotated_transcripts": panel["transcript_count"],
            "exact_cdna_classes": panel["equivalence_group_count"],
            "selected_assays": panel["selected_assay_count"],
            "uncovered_classes": len(panel["uncovered_equivalence_group_ids"]),
            "unresolved_end_pairs": len(panel["unresolved_group_pairs"]),
            "genomic_specificity": (
                "completed"
                if panel.get("genomic_specificity_assessments")
                else "not_run"
            ),
        }
    )
    for offset, assay in enumerate(panel["selected_assays"], start=1):
        label = f"A{offset}"
        pair = assay["primer_pair"]
        test = read_json(ROOT / f"{gene}_{label}_test.json")
        primer_rows.append(
            {
                "gene": gene,
                "assay": label,
                "source": "GENtle/Primer3 end-matrix design",
                "forward_5to3": pair["forward"]["sequence"],
                "reverse_5to3": pair["reverse"]["sequence"],
                "forward_tm_c": f"{pair['forward']['tm_c']:.2f}",
                "reverse_tm_c": f"{pair['reverse']['tm_c']:.2f}",
                "design_transcript": assay["design_transcript_id"],
                "assay_id": assay["assay_id"],
                "oligo_qc": test["oligo_qc"]["status"],
                "genomic_specificity": "not_run",
            }
        )
        for transcript in test["transcript_results"]:
            for product in transcript["products"]:
                f_hit = transcript["forward_hits"][product["forward_hit_index"]]
                r_hit = transcript["reverse_hits"][product["reverse_hit_index"]]
                binding_rows.append(
                    {
                        "gene": gene,
                        "assay": label,
                        "transcript_id": transcript["transcript_id"],
                        "transcript_label": transcript["transcript_label"],
                        "cdna_length_bp": transcript["cdna_length_bp"],
                        "forward_cdna_1based": f"{f_hit['start_0based'] + 1}-{f_hit['end_0based_exclusive']}",
                        "forward_exon": exon_label(f_hit, transcript["exon_segments"]),
                        "reverse_cdna_1based": f"{r_hit['start_0based'] + 1}-{r_hit['end_0based_exclusive']}",
                        "reverse_exon": exon_label(r_hit, transcript["exon_segments"]),
                        "locus_local_forward_1based": source_ranges(f_hit),
                        "locus_local_reverse_1based": source_ranges(r_hit),
                        "product_bp": product["amplicon_length_bp"],
                        "spans_junction": product["spans_junction"],
                        "genomic_carryover_risk": product["genomic_carryover_risk"],
                    }
                )

legacy = (
    (
        "L1",
        "existing e2f/e3r",
        "TCAAGCAGGTGCACACTTCTGAG",
        "TGGGACGACCTCCACAAAGC",
        "PATZ1_existing_e2f_e3r_test.json",
    ),
    (
        "L2",
        "existing e4f/e6r",
        "GGCCCAGCAACTTCTGCAGTATC",
        "TGAGAGGTCACCATAGGAGTCAGAG",
        "PATZ1_existing_e4f_e6r_test.json",
    ),
)
for label, source, forward, reverse, filename in legacy:
    test = read_json(ROOT / filename)
    primer_rows.append(
        {
            "gene": "PATZ1",
            "assay": label,
            "source": source,
            "forward_5to3": forward,
            "reverse_5to3": reverse,
            "forward_tm_c": "",
            "reverse_tm_c": "",
            "design_transcript": "",
            "assay_id": "",
            "oligo_qc": test["oligo_qc"]["status"],
            "genomic_specificity": "not_run",
        }
    )
    for transcript in test["transcript_results"]:
        for product in transcript["products"]:
            f_hit = transcript["forward_hits"][product["forward_hit_index"]]
            r_hit = transcript["reverse_hits"][product["reverse_hit_index"]]
            binding_rows.append(
                {
                    "gene": "PATZ1",
                    "assay": label,
                    "transcript_id": transcript["transcript_id"],
                    "transcript_label": transcript["transcript_label"],
                    "cdna_length_bp": transcript["cdna_length_bp"],
                    "forward_cdna_1based": f"{f_hit['start_0based'] + 1}-{f_hit['end_0based_exclusive']}",
                    "forward_exon": exon_label(f_hit, transcript["exon_segments"]),
                    "reverse_cdna_1based": f"{r_hit['start_0based'] + 1}-{r_hit['end_0based_exclusive']}",
                    "reverse_exon": exon_label(r_hit, transcript["exon_segments"]),
                    "locus_local_forward_1based": source_ranges(f_hit),
                    "locus_local_reverse_1based": source_ranges(r_hit),
                    "product_bp": product["amplicon_length_bp"],
                    "spans_junction": product["spans_junction"],
                    "genomic_carryover_risk": product["genomic_carryover_risk"],
                }
            )


def write_tsv(name: str, rows: list[dict]) -> None:
    with (ROOT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


write_tsv("gene_summary.tsv", summary_rows)
write_tsv("primer_pairs.tsv", primer_rows)
write_tsv("isoform_binding_positions.tsv", binding_rows)

# A compact isoform-by-assay view complements the long coordinate table and
# mirrors the virtual gel: assays are columns, while each non-empty cell gives
# the exact predicted product and transcript-local primer coordinates.
matrix_rows: list[dict] = []
for gene in GENES:
    gene_primers = [row["assay"] for row in primer_rows if row["gene"] == gene]
    transcripts = sorted(
        {
            (row["transcript_id"], row["transcript_label"])
            for row in binding_rows
            if row["gene"] == gene
        }
    )
    for transcript_id, transcript_label in transcripts:
        matrix_row = {
            "gene": gene,
            "transcript_id": transcript_id,
            "transcript_label": transcript_label,
        }
        for assay in gene_primers:
            products = [
                row
                for row in binding_rows
                if row["gene"] == gene
                and row["assay"] == assay
                and row["transcript_id"] == transcript_id
            ]
            matrix_row[assay] = "; ".join(
                f"{row['product_bp']} bp (F {row['forward_cdna_1based']}; "
                f"R {row['reverse_cdna_1based']})"
                for row in products
            )
        matrix_rows.append(matrix_row)

with (ROOT / "isoform_assay_matrix.tsv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    fieldnames = ["gene", "transcript_id", "transcript_label"] + [
        f"{gene}:{row['assay']}"
        for gene in GENES
        for row in primer_rows
        if row["gene"] == gene
    ]
    # Use a common rectangular schema so all five genes remain in one file.
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in matrix_rows:
        expanded = {
            "gene": row["gene"],
            "transcript_id": row["transcript_id"],
            "transcript_label": row["transcript_label"],
        }
        for assay, value in row.items():
            if assay not in {"gene", "transcript_id", "transcript_label"}:
                expanded[f"{row['gene']}:{assay}"] = value
        writer.writerow(expanded)

for gene in GENES:
    assays = [row["assay"] for row in primer_rows if row["gene"] == gene]
    rows = [row for row in matrix_rows if row["gene"] == gene]
    with (ROOT / f"{gene}_isoform_assay_matrix.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["transcript_id", "transcript_label", *assays],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def md_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    lines.extend("| " + " | ".join(str(cell) for cell in row) + " |" for row in rows)
    return "\n".join(lines)


primer_table = md_table(
    ["Gene", "Lane", "Forward primer (5′→3′)", "Reverse primer (5′→3′)", "QC", "Origin"],
    [
        [row["gene"], row["assay"], f"`{row['forward_5to3']}`", f"`{row['reverse_5to3']}`", row["oligo_qc"], row["source"]]
        for row in primer_rows
    ],
)
summary_table = md_table(
    ["Gene", "Ensembl transcripts", "Exact cDNA classes", "Assays", "Panel", "Uncovered classes"],
    [
        [
            row["gene"],
            row["annotated_transcripts"],
            row["exact_cdna_classes"],
            row["selected_assays"],
            row["completion_status"],
            row["uncovered_classes"],
        ]
        for row in summary_rows
    ],
)

readme = f"""# TP73-cofactor transcript primer panels

This is an in-silico design report, not yet an ordering recommendation. GENtle
used live Ensembl release 116 / GRCh38 gene models, Primer3 2.6.1, exact cDNA
matching, a 10 kb endpoint-product ceiling, and oligo-dT cDNA as the experimental
context. The four genes other than PATZ1 have Clariom D prioritisation but no
independent transcript evidence in the current project; predicted compatibility
with an annotated transcript is therefore not evidence that the transcript is
expressed in the cells.

## Coverage summary

{summary_table}

AHCTF1 and POU2F2 are deliberately marked partial. A partial panel must not be
described as complete merely because the designed lanes are informative.
AHCTF1 covers all 18 exact cDNA classes but does not realize every requested
end-class reaction. POU2F2 leaves the short POU2F2-206 (785 nt) and POU2F2-209
(700 nt) cDNA classes uncovered.

## Primer pairs and gel lanes

{primer_table}

For PATZ1, lanes L1 and L2 are the existing e2f/e3r and e4f/e6r assays. The six
GENtle-designed lanes are A1–A6. For the other genes, lanes are A1 upward.

## Exact transcript binding

`isoform_binding_positions.tsv` gives one row per predicted product. cDNA primer
positions are 1-based inclusive coordinates in the named Ensembl transcript;
exon ordinals are transcript-local. Locus-local coordinates are 1-based within
the imported Ensembl gene interval, not chromosome coordinates. The complete
GENtle reports retain exon segments and the imported locus provenance.
`isoform_assay_matrix.tsv` presents the same results as an isoform-by-primer-pair
matrix; each non-empty cell contains product size and forward/reverse cDNA
coordinates. The five `{{GENE}}_isoform_assay_matrix.tsv` files provide the same
view with only that gene's assay columns.

## GENtle virtual gels

- `PATZ1_virtual_gel_complete.svg` — A1–A6 plus existing L1–L2
- `TFAP2C_virtual_gel.svg`
- `AHCTF1_virtual_gel.svg`
- `POU2F2_virtual_gel.svg`
- `LMX1B_virtual_gel.svg`

### PATZ1

![PATZ1 GENtle virtual gel](PATZ1_virtual_gel_complete.svg)

### TFAP2C

![TFAP2C GENtle virtual gel](TFAP2C_virtual_gel.svg)

### AHCTF1

![AHCTF1 GENtle virtual gel](AHCTF1_virtual_gel.svg)

### POU2F2

![POU2F2 GENtle virtual gel](POU2F2_virtual_gel.svg)

### LMX1B

![LMX1B GENtle virtual gel](LMX1B_virtual_gel.svg)

Each sample lane is one primer pair. GENtle materialized every predicted cDNA
product as a sequence and used its shared serial gel renderer. Close products
that would co-migrate are drawn as one apparent band but remain separate in
`isoform_binding_positions.tsv` and the assay JSON.

GENtle is therefore capable of producing the requested figures. Its public
materialization operation currently repeats the gene name as every lane label,
so the displayed A1/A2/... and PATZ1 L1/L2 labels were applied in a
visualization-only copy of the project state. No sequence, predicted product,
arrangement, or gel parameter was changed.

## Important interpretation boundaries

- PATZ1 protein-size evidence does not by itself identify one unique Ensembl
  transcript. L1 produces the same 142 bp product from PATZ1-201 and PATZ1-205;
  L2 detects several other annotated transcripts with products near 92–233 bp.
- Oligo-dT reverse transcription may underrepresent distant 5′ sequence. This is
  especially serious for AHCTF1; failure of a long product cannot establish
  transcript absence.
- Clariom D probe/junction behavior is prioritisation evidence, not transcript
  quantification.
- Five proposed lanes failed the conservative exact-run oligo screen and should
  be redesigned before ordering: PATZ1 A1/A3 (5-base homopolymer), AHCTF1 A1
  (6-base 3′ self-complementarity), POU2F2 A8 (5-base homopolymer), and POU2F2
  A14 (4-base 3′ self-complementarity). The two legacy PATZ1 pairs passed this
  screen.
- Whole-genome specificity is reported separately. Until it passes, these
  sequences must not be treated as order-ready.
- A whole-GRCh38 preparation was started but stopped before indexing: the
  expanded Ensembl FASTA and GTF consumed about 7.5 GB and left the host at 99%
  filesystem use. Consequently genomic specificity remains `not_run`; it was
  not converted into a vacuous pass.
"""
(ROOT / "README.md").write_text(readme, encoding="utf-8")

# GENtle's current materialization operation names every lane container after
# the gene.  Make a visualization-only project copy with concise lane names;
# sequences, products, arrangements, and all biological results are unchanged.
state_path = ROOT / "five_genes.with_products.json"
if state_path.exists():
    state = read_json(state_path)
    for container in state["container_state"]["containers"].values():
        members = container.get("members", [])
        if not members:
            continue
        first = members[0]
        if first.startswith("PATZ1_existing_e2f_e3r_"):
            container["name"] = "L1 e2f/e3r"
        elif first.startswith("PATZ1_existing_e4f_e6r_"):
            container["name"] = "L2 e4f/e6r"
        else:
            for gene in GENES:
                prefix = f"{gene}_A"
                if first.startswith(prefix):
                    lane = first[len(gene) + 1 :].split("_", 1)[0]
                    container["name"] = lane
                    break
    with (ROOT / "five_genes.gel_labels.state.json").open("w", encoding="utf-8") as handle:
        json.dump(state, handle, separators=(",", ":"))
