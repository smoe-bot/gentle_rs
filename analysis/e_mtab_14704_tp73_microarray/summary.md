# E-MTAB-14704 Targeted First-Pass Summary

Dataset: SK-MEL-29 Clariom D arrays, 3 replicates each for Ad_DNp73beta,
Ad_GFP, and Ad_TAp73alpha.

Processing: APT `apt-cel-extract` feature intensities, Clariom D Bioconductor
platform annotation, quantile normalization across the nine arrays, log2 values,
median summarization over probes. Statistics are exploratory Welch tests with
Benjamini-Hochberg correction across the targeted set. Paired t-test p/q columns
are also reported, assuming replicate suffixes 1/2/3 are matched experimental
rounds; the SDRF sample sheet does not explicitly document that block.

## TP73

The TP73 core transcript cluster is strongly elevated in both TP73 constructs
relative to GFP:

- GFP intrinsic TP73: mean log2 6.55
- DNp73beta: mean log2 12.8, log2FC vs GFP 6.24, q 0.000171
- TAp73alpha: mean log2 13.3, log2FC vs GFP 6.70, q 0.00966
- TAp73alpha vs DNp73beta: log2FC 0.46, q 0.229

With the replicate-suffix paired assumption, the TP73 construct effects remain
strong:

- TAp73alpha vs GFP: paired q 0.0187
- DNp73beta vs GFP: paired q 0.00288

The regional probe/probeset pattern is more informative than the collapsed core
cluster. Several TP73 probe-set regions are high in both constructs, while
others are strongly TAp73alpha-specific and remain near the GFP baseline in
DNp73beta. The strongest TAp73alpha-specific regions include:

- PSR0100145789.hg.1: TAp73alpha +7.82 log2 vs GFP; DNp73beta +0.19
- PSR0100145802.hg.1: TAp73alpha +7.09; DNp73beta +0.01
- PSR0100145782.hg.1: TAp73alpha +6.80; DNp73beta -0.12
- PSR0100145783.hg.1: TAp73alpha +6.73; DNp73beta -0.11
- PSR0100145805.hg.1: TAp73alpha +5.91; DNp73beta -0.06

Regions high in both constructs include PSR0100145793.hg.1,
PSR0100145792.hg.1, PSR0100145794.hg.1, PSR0100145801.hg.1,
PSR0100145800.hg.1, PSR0100145796.hg.1, and PSR0100145797.hg.1.

## Lists A and B

Among List A/B genes, the clearest FDR-supported TAp73alpha-vs-GFP changes are:

- TFAP2C: +2.46 log2, q 0.0295
- ELK3: +1.20 log2, q 0.0295

Nominal but not FDR-supported List A/B changes include:

- TFAP2A: -0.85 log2, p 0.021, q 0.133
- INSM1: +0.85 log2, p 0.029, q 0.134
- FOXF2: +0.24 log2, p 0.025, q 0.133
- GLIS3: +0.18 log2, p 0.043, q 0.153

No List A/B gene is FDR-supported for DNp73beta-vs-GFP in this targeted
first-pass analysis.

Under the paired-suffix sensitivity analysis, List A/B support weakens rather
than strengthens after correction:

- TFAP2C: paired p 0.0053, paired q 0.085
- ELK3: paired p 0.0166, paired q 0.174
- INSM1: paired p 0.0218, paired q 0.174

So the robust interpretation is that TFAP2C and ELK3 remain the main List A
signals by effect size and nominal/uncorrected support, while formal FDR support
depends on using the unpaired Welch model.

## TP73-FUS Panel

No TP73-FUS-panel gene reaches FDR support in this first-pass targeted analysis.
Notable nominal or directional observations:

- HDAC1: TAp73alpha vs GFP +0.88 log2, p 0.039, q 0.153
- HDAC2: TAp73alpha vs GFP -0.42 log2, p 0.020, q 0.133
- MDM2: DNp73beta vs GFP -0.88 log2, TAp73alpha vs GFP +0.83, and
  TAp73alpha vs DNp73beta +1.71; the TA-vs-DN nominal p is 0.0059 but q 0.11
- FUS: essentially unchanged vs GFP; TAp73alpha is slightly lower than
  DNp73beta (-0.22 log2), nominal p 0.043, q 0.198
- E2F1: modest DNp73beta elevation (+0.55 log2 vs GFP), not significant

The paired-suffix sensitivity analysis does not add FDR-supported TP73-FUS panel
hits. The most interesting nominal paired patterns remain HDAC1, HDAC2, MDM2,
REST, and TARDBP, but all paired q-values are above 0.18 in the targeted set.

Output tables:

- `target_gene_summary.csv`
- `list_a_b_gene_summary.csv`
- `tp73_fus_panel_gene_summary.csv`
- `tp73_regional_probe_pattern.csv`
- `tp73_top_changing_probe_regions.csv`

## Oligo/Limma Second Pass

After R, `oligo`, `limma`, and a local `pd.clariom.d.human` library became
available, a whole-array RMA/limma pass was added:

- `oligo_rma_core_expression.csv`
- `oligo_limma_unpaired_all_core.csv`
- `oligo_limma_paired_all_core.csv`
- `oligo_limma_unpaired_target_transcripts.csv`
- `oligo_limma_paired_target_transcripts.csv`

These limma adjusted p-values are across the whole core transcript-cluster
array, not only across the small target list.

A targeted-panel FDR table is also written:

- `oligo_limma_paired_targeted_panel_fdr.csv`

This uses paired limma p-values but applies Benjamini-Hochberg correction only
inside the pre-specified panels.

For TAp73alpha vs GFP, the paired limma model supports:

- TP73: logFC +7.14, adjusted p 1.03e-05
- TFAP2C: logFC +3.01, adjusted p 0.000115
- INSM1: logFC +1.83, adjusted p 0.00341
- EGR3: logFC +1.38, adjusted p 0.00650
- ELK3: logFC +1.41, adjusted p 0.0105
- TFAP2A: logFC -0.80, adjusted p 0.0469

For DNp73beta vs GFP, only TP73 is supported after whole-array correction in
the paired limma model:

- TP73: logFC +6.84, adjusted p 0.000418

For TAp73alpha vs DNp73beta, the paired limma model supports:

- TFAP2C: logFC +2.28, adjusted p 0.000594
- INSM1: logFC +1.70, adjusted p 0.00455
- MDM2: logFC +1.60, adjusted p 0.00683
- EGR3: logFC +1.10, adjusted p 0.0185
- HDAC2: logFC -0.74, adjusted p 0.0294

This makes the biological summary cleaner: DNp73beta primarily shows the TP73
construct signal in the target set, while TAp73alpha has a stronger downstream
transcriptional signature involving TFAP2C, INSM1, EGR3, ELK3, TFAP2A, MDM2,
and HDAC2.

For a targeted FUS/ALS-associated panel FDR, the interpretation is less strict:

- TAp73alpha vs GFP: HDAC1 (+0.79, targeted q 0.010), HDAC2 (-0.62,
  targeted q 0.010), and MDM2 (+0.70, targeted q 0.046) pass panel-level FDR.
- DNp73beta vs GFP: MDM2 (-0.90, targeted q 0.037) passes panel-level FDR.
- TAp73alpha vs DNp73beta: MDM2 (+1.60, targeted q 0.00073), HDAC2 (-0.74,
  targeted q 0.0034), and HDAC1 (+0.55, targeted q 0.040) pass panel-level FDR.

PATZ1 is present in List A but is not an upregulated expression hit. In the
paired limma model it trends downward for TAp73alpha:

- TAp73alpha vs GFP: logFC -0.62, whole-array adjusted p 0.436
- TAp73alpha vs DNp73beta: logFC -0.85, targeted List A q 0.050 and
  whole-array adjusted p 0.205
