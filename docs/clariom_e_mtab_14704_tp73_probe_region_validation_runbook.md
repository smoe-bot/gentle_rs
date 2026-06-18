# E-MTAB-14704 TP73 Probe-Region Validation Runbook

This runbook documents the explicit local workflow for validating the Clariom D
probe-region bridge on the public E-MTAB-14704 TP73 context. GENtle must not
silently download CEL files, install R/Bioconductor packages, run APT, or fetch
vendor annotation payloads.

The committed validation fixture at
`test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation/`
is generated from a compact Glen-style adapter input with TP73 GRCh38.p14
coordinates. It exists so CI, GUI tests, and release review can exercise
projection and interpretation surfaces without raw arrays or Glen's full
exploratory analysis snapshot.

## Inputs

- Publication-resource dataset:
  `data/publication_resources/rostock_p73_clariomd_e_mtab_14704/`
  - `manifest.json`
  - `download_manifest.tsv`
  - `download.sh`
  - `E-MTAB-14704.idf.txt`
  - `E-MTAB-14704.sdrf.txt`
- TP73 anchor: `test_files/tp73.ncbi.gb`
  - GRCh38.p14
  - chromosome 1, 3652516..3736201
- Optional real execution tools, installed and run manually by the developer:
  - Rscript
  - Bioconductor `oligo`, `limma`, and `pd.clariom.d.human`
  - or Affymetrix Power Tools plus compatible Clariom D library files
- Optional vendor support ZIPs:
  `data/resources/affymetrix/clariom_d_human_na36_hg38/`
- Compact adapter input:
  `test_files/fixtures/probe_region_adapter_inputs/clariom_e_mtab_14704_tp73_glen_style.csv`

## Preflight The Dataset

Inspect the committed resource metadata:

```bash
cargo run --bin gentle_cli -- resources status-publication-dataset E-MTAB-14704
```

Plan the TP73 probe-region analysis without running external tools:

```bash
cargo run --bin gentle_cli -- arrays probe-regions \
  --dataset E-MTAB-14704 \
  --gene TP73 \
  --platform Clariom_D_Human \
  --normalization rma \
  --output analysis/probe_regions/e_mtab_14704_tp73 \
  --dry-run
```

Expected behavior:

- the plan resolves the nine declared local CEL paths from the resource catalog
- locally present SDRF metadata is used with sample column `Array Data File`
  and condition column `Characteristics[genetic modification]`
- missing CEL files, R packages, APT binaries, or vendor support files are
  reported as preflight findings
- `analysis/probe_regions/e_mtab_14704_tp73/plan.json` is written as the
  inspectable versioned plan artifact
- no bytes are downloaded and no analysis backend is executed

## Explicit Local Raw-Data Preparation

Only run this when the local machine should fetch the raw CEL files:

```bash
cd data/publication_resources/rostock_p73_clariomd_e_mtab_14704
./download.sh
```

Alternatively, use the resource helper explicitly:

```bash
cargo run --bin gentle_cli -- resources prepare-publication-dataset \
  E-MTAB-14704 \
  --categories raw_microarray \
  --download-files
```

Raw CEL files remain uncommitted.

## Manual Backend Execution

Use the generated preflight command as the starting point for a local R/oligo
run. The generic helper is:

```bash
Rscript scripts/probe_regions_oligo.R \
  --cel data/publication_resources/rostock_p73_clariomd_e_mtab_14704/P_SKMel29_AdGFP_1.CEL \
  --metadata data/publication_resources/rostock_p73_clariomd_e_mtab_14704/E-MTAB-14704.sdrf.txt \
  --sample-column "Array Data File" \
  --condition-column "Characteristics[genetic modification]" \
  --gene TP73 \
  --platform-package pd.clariom.d.human \
  --normalization rma \
  --output analysis/probe_regions/e_mtab_14704_tp73
```

Add all nine CEL files for the real analysis. If using APT instead, run APT
manually and then import its summary, annotation, metadata, and PM-probe
intensity tables with:

```bash
cargo run --bin gentle_cli -- arrays import-apt-probe-region-output \
  SUMMARY.tsv \
  ANNOTATION.csv \
  analysis/probe_regions/e_mtab_14704_tp73 \
  --metadata data/publication_resources/rostock_p73_clariomd_e_mtab_14704/E-MTAB-14704.sdrf.txt \
  --condition-column "Characteristics[genetic modification]" \
  --sample-column "Array Data File" \
  --probe-intensity PROBE_INTENSITY.tsv \
  --probe-id-column probe_id \
  --platform Clariom_D_Human \
  --normalization rma \
  --coordinate-system hg38 \
  --genome-build GRCh38.p14
```

## Probe-Set Activity Presentation

For figure preparation from a local full E-MTAB-14704/Affymetrix working set,
GENtle keeps the code path in Git while leaving bulky CEL-derived
intermediates outside the repository. The helper below renders a selected
gene-panel view from explicit local inputs:

```bash
python3 scripts/render_clariomd_probe_set_activity.py \
  --raw-features analysis/e_mtab_14704_tp73_microarray/all_arrays_raw_features.tsv \
  --sqlite analysis/e_mtab_14704_tp73_microarray/Rlib/pd.clariom.d.human/extdata/pd.clariom.d.human.sqlite \
  --vendor-probeset-zip data/publication_resources/rostock_p73_clariomd_e_mtab_14704/library/Clariom_D_Human-na36-hg38-probeset-csv.zip \
  --output-dir analysis/e_mtab_14704_tp73_microarray/gene_panel_probe_set_activity \
  --genes TP73,FUS,PATZ1,E2F1,TARDBP,PLK1,TERT,HDAC1,HDAC2,HDAC6
```

Expected local-only outputs include:

- `probe_set_activity_summary.tsv`: probeset-level raw PM-probe means and
  group-level log2 contrasts.
- `probe_level_activity.tsv`: selected PM-probe raw intensities.
- `gene_contrast_probe_set_summary.png/.svg/.pdf`: compact per-gene contrast
  distributions.
- `probe_set_individual_arrays_heatmap_10_gene.png/.pdf`: the nine arrays as
  individual columns, grouped by paired experiment/time/person.
- `probe_set_paired_contrast_heatmap_10_gene.png/.pdf`: within-experiment
  contrasts such as `TAp73alpha_i - GFP_i` and `DNp73beta_i - GFP_i`.
- `paired_gene_level_summary.tsv`: per-gene median paired contrasts.

These outputs are deliberately uncommitted derived analysis artifacts. The raw
feature table can be regenerated from CEL files with APT-style probe extraction
or equivalent local tooling; the script only consumes that explicit table and
does not claim a formal expression model or isoform-support verdict.

## Committed Fixture Validation Path

The committed validation fixture is deterministic and safe for CI:

```bash
cargo run --bin gentle_cli -- arrays inspect-probe-region-output \
  test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation
```

The fixture is regenerated by
`GentleEngine::import_glen_probe_region_fixture` from the compact adapter input.
The committed output has three parent region rows and fourteen PM-probe rows,
all marked `probe_level_input`. The selected rows exercise shared exon overlap,
unique exon overlap, one junction-spanning interval, and one non-exonic
transcript-span constraint.

Project the PM-probe rows into a TP73-anchored state or workflow. The example
below assumes the active GENtle state already contains an `array_slice` sequence
with TP73-compatible genome-anchor metadata:

```bash
cargo run --bin gentle_cli -- arrays project-probe-region-output \
  array_slice \
  test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation \
  --level pm_probe \
  --contrast AdTAp73alpha-AdGFP \
  --min-abs-logfc 0.5 \
  --max-features 20 \
  --clear-existing
```

Interpret the projected rows as a review-only transcript/exon geometry report:

```bash
cargo run --bin gentle_cli -- arrays interpret-probe-region-evidence \
  array_slice \
  --gene TP73 \
  --level pm_probe \
  --min-abs-logfc 0.5 \
  --path artifacts/clariom_e_mtab_14704_tp73_interpretation.json
```

The interpretation report is a constraint/triage artifact only. It records
geometry compatibility, shared-vs-unique overlap, junction spans, and ambiguity
tags such as `probe_sequence_alignment_not_assessed`,
`multi_hit_not_assessed`, and `isoform_support_not_inferred`. It is not a probe
specificity verdict, multi-hit assessment, or isoform-support call.

## Deterministic Tests

Run the focused validation checks:

```bash
cargo test --lib -- \
  glen_probe_region_adapter_regenerates_committed_e_mtab_fixture \
  main_area_dna_projects_probe_region_output_through_shared_shell_capability \
  main_area_dna_interprets_probe_region_evidence_through_shared_shell_capability \
  e_mtab_14704_tp73_validation_fixture_projects_and_interprets_pm_probe_evidence \
  e_mtab_14704_tp73_validation_report_is_probe_location_figure_ready
```

Useful neighboring checks:

```bash
cargo test --lib -- interpret_probe_region_evidence
cargo test --lib -- project_probe_region_output_pm_probe
cargo test --lib -- import_apt_probe_region_output
```

Before handoff:

```bash
cargo check -q
git diff --check
```
