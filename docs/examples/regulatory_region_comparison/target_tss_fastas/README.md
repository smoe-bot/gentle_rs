# Target-gene TSS FASTA export

This bundle exports one FASTA file per target gene from a validated GENtle
transcript promoterome. Each record represents one distinct genomic TSS and
strand. Transcripts sharing that exact TSS are retained together in the header.

The requested sequence is **−500 through +200 bp in transcript orientation**,
including the TSS (701 bp). Thus every record is written 5′→3′ relative to its
transcript, including reverse-complemented sequence for minus-strand genes.
Coordinates in headers remain GRCh38 1-based inclusive genomic coordinates.

## Replay

Run from the repository root at the producer revision recorded in
`manifest.json`:

```bash
python3 scripts/export_target_tss_fastas.py \
  --promoterome /path/to/verified/human-grch38-ensembl116-promoterome \
  --expected-genome-id 'Human GRCh38 Ensembl 116' \
  --source-revision "$(git rev-parse HEAD)" \
  --gene CD44 --gene TGFB1 --gene SERPINE1 --gene PATZ1 --gene TP73 \
  --output docs/examples/regulatory_region_comparison/target_tss_fastas/bundle
```

The exporter validates the promoterome receipt and all consumed source hashes,
checks unclipped strand/TSS geometry, slices the already transcript-oriented
source sequence, and records source and output digests in `manifest.json` and
`SHA256SUMS`.

Transcript annotation does not establish which TSS is used in the assayed
cells, and a TSS window is not proof of promoter activity.
