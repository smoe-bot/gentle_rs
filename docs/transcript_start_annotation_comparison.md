# Transcript-start annotation comparison

GENtle can prepare the current catalogued human NCBI RefSeq assembly without a
separate download helper:

```sh
gentle_cli genomes prepare \
  "Human GRCh38 NCBI RefSeq GCF_000001405.40" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes \
  --timeout-secs 7200

gentle_cli genomes status \
  "Human GRCh38 NCBI RefSeq GCF_000001405.40" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

The catalog derives the official NCBI genomic FASTA and genomic GFF3 URLs from
assembly accession `GCF_000001405.40` and assembly name `GRCh38.p14`. Preparation
downloads both files, records their source identities and hashes, and builds the
FASTA, gene, transcript and BLAST indexes. `GCF_` denotes the RefSeq assembly;
the paired GenBank assembly is `GCA_000001405.29`. They are not interchangeable
identifiers, and the report does not relabel RefSeq as GenBank.

## Compare with a bound Ensembl TSS report

`scripts/compare_transcript_start_annotations.py` compares one existing
`gentle.tss_tfbs_profiles.v1` report with the transcript index from a prepared
secondary annotation:

```sh
python3 scripts/compare_transcript_start_annotations.py \
  --primary-tss-report /path/to/report.json \
  --secondary-manifest data/genomes/human_grch38_ncbi_refseq_gcf_000001405_40/manifest.json \
  --catalog assets/genomes.json \
  --gene CD44 --gene TGFB1 --gene SERPINE1 --gene PATZ1 --gene TP73 \
  --producer-revision "$(git rev-parse HEAD)" \
  --gentle-cli /path/to/exact/gentle_cli \
  --output-report /fresh/output/annotation-comparison.json \
  --output-dir /fresh/output/pages \
  --output-receipt /fresh/output/annotation-comparison.receipt.json
```

The result is `gentle.transcript_start_annotation_comparison.v1`. For each gene
it retains every source TSS coordinate and transcript membership, exact-coordinate
agreement counts, nearest opposite-source coordinates with signed
transcript-direction deltas, and the original reporter-selection flag.

Positive delta means that the RefSeq start lies downstream from the Ensembl
start in transcript direction; negative means upstream. The sign is therefore
reversed on a minus-strand gene relative to increasing genomic coordinates.
Nearest-coordinate relationships are descriptive. They do not establish
transcript orthology, equivalence, a consensus TSS, or TSS use in the assayed
cells.

The comparison receipt binds the primary TSS report, catalog, prepared-genome
manifest, NCBI annotation, derived transcript index, canonical JSON report, and
one aligned SVG page per gene. The NCBI annotation SHA-1 is rechecked against
GENtle's install manifest and an additional SHA-256 is recorded. Missing genes,
mixed strands/chromosomes, assembly-core disagreement and changed files fail
closed.

## Reporter-report integration

The integrated locus/TSS compositor accepts an optional comparison triple:

```text
--annotation-comparison-report  annotation-comparison.json
--annotation-comparison-receipt annotation-comparison.receipt.json
--annotation-comparison-svg     pages/CD44_Ensembl_RefSeq_TSS_comparison.svg
```

When supplied together, the comparison becomes page 2, between the locus
overview and selected-TSS TFBS pages. Its 1400-pixel page and x=255..1050 plot
frame match the other pages. The compositor verifies the comparison's gene,
primary reference, primary-report hash, page hash and geometry before rendering.
The original Ensembl models, CUT&RUN context, TFBS arrays and reporter selection
remain unchanged.

This is the first consumer of a source-neutral annotation-comparison artifact.
Other reports may reuse the same JSON, but should not silently merge annotations
or promote proximity to biological agreement.
