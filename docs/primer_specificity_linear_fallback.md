# Linear fallback for exhaustive primer-specificity handoffs

`scripts/primer_specificity_linear_fallback.mjs` is an experimental,
fail-closed fallback for a specific scalability boundary: an exhaustive BLAST
handoff can contain many thousands of HSP rows per primer, while the current
GENtle importer fetches a subject window and performs semiglobal realignment
for every row before grouping candidate products.

The fallback consumes existing `gentle.primer_specificity_handoff.v1` or
`gentle.transcript_assay_panel_specificity_handoff.v1` files. It reruns only
their declared `blastn` commands, changing the output path and adding
`qlen qseq sseq` to the tabular output. It then:

1. derives full-query coverage, terminal mismatches, and total mismatches from
   the aligned strings;
2. filters hits with the handoff's recorded coverage and 3-prime policy;
3. groups hits by subject before pairing inward-facing primers;
4. applies the handoff's readiness amplicon ceiling;
5. requires every declared intended cDNA subject and rejects undeclared
   products below the recorded mismatch threshold.

This route is deliberately labelled a fallback. Its result is independently
auditable BLAST-on-cDNA evidence, not a `GENtle` finalized specificity pass.
It exists so large evidence sets can be studied and compared while the native
importer is optimized.

## Request

```json
{
  "schema": "gentle.primer_specificity_linear_fallback_request.v1",
  "handoff_paths": [
    "/absolute/path/panel_specificity_handoff.json",
    "/absolute/path/primer_specificity_handoff.json"
  ],
  "output_dir": "/new/absolute/output/directory"
}
```

The output directory must not already exist. Run:

```sh
node scripts/primer_specificity_linear_fallback.mjs request.json
```

Tests use only Node's standard library:

```sh
node --test scripts/primer_specificity_linear_fallback.test.mjs
```

## Integration direction

The preferred long-term fix is engine-owned: retain exhaustive database search
completeness while avoiding one indexed-sequence fetch and dynamic-programming
alignment per raw HSP. This prototype makes the aligned-string and
subject-grouping behavior explicit enough to compare against that native
implementation.
