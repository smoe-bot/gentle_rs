# Bounded Post-BLAST Primer-Import Phase Attribution

Executed on 2026-11-05 following `docs/glen_post_blast_import_phase_trace.md`.
The tested GENtle source was the clean detached revision
`dfca32611ad2b217382186038a4c3879d907fdde`, tree
`b5428910481eeaf23fe23ec0e04d885997fe62ac`. Its parent was
`108cc5e51d02b31f7af90d941821a2733ad50059`; their diff changed only four
documentation files. No BLAST search, reference preparation, source repair,
profile override, profiler, package, GUI, or biological acceptance was run.

## Source, host, and binary

- `Cargo.lock` SHA-256:
  `7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.
- Linux x86_64, kernel `7.2.7+deb14-amd64`; Rust/Cargo 1.97.1.
- Ordinary dev/default-feature CLI, one Cargo job, `CARGO_INCREMENTAL=0`.
  No opt-level, frame-pointer, `gui-profiler`, or package-profile override.
- CLI: 878,359,560 bytes; SHA-256
  `b07e5b0d941bf6e9576cfa2789ba36e220b6f43b01667d664aba0115f5227ec4`.
- `--version` identifies source revision
  `0.1.0-internal.13+git.dfca32611ad2b217382186038a4c3879d907fdde`.

## Admission

- Original audit archive: 103,172,261 bytes; SHA-256
  `35e000d086f17348ea76dbaf88712f44b808abea045a72509cabc6d6034a1dcd`.
- All 213 archived manifest hashes passed before and after the attempts.
- ACTB: 5,999 HSP rows / 686,525 bytes. TP53: 12,202 rows / 1,404,773
  bytes. Each case retained four completed, exit-zero command outputs whose
  byte sizes and SHA-256 hashes matched its manifest.
- The moved Ensembl-116 cDNA BLAST index was restored byte-identically from
  the verified Storage Box archive to its content-bound original path. The ten
  restored files passed the archive's SHA-256 manifest. No alias, JSON rewrite,
  or database preparation was used.
- Exact-source pre- and post-run inspections were byte-identical:
  - cDNA fingerprint
    `sha256:3159519891bd8c0d4690e9c2f3de8eeb9621152b58c749e3d7d4920f4913263f`,
    465,769 sequences / 1,276,790,359 bases; annotation SHA-1
    `5f571d5775368aeb4e8a88c3ce4b1c8c138e3acf`; unmasked transcriptome index.
  - GRCh38 fingerprint
    `sha256:90a64ae8ff24fdbab67286bab563c95a6b2a04be37f50ad811474106c1ba3e7e`,
    706 contigs / 3,291,585,349 bases; genomic index.

## Focused checks

All requested commands passed, with complete stdout/stderr/time records:

- `primer_import_trace`: 5 passed, 0 failed, 0 ignored; 4,941 filtered out.
- `engine::primer_specificity_multi::tests::`: 10 passed, 0 failed,
  0 ignored; 4,936 filtered out.
- `specificity_matrix_`: 13 passed, 0 failed, 0 ignored; 4,933 filtered out.
- `cargo check -q --locked`: passed.

The first cold test build took 20:42.91 wall and 12,386,816 KiB peak RSS.
The exact CLI build took 7:51.37 wall and 10,129,224 KiB peak RSS.

## Bounded real attempts

Each attempt used a distinct initially absent state. Commands and trace
environment are retained per attempt. No attempt left a state, stdout,
run-owned process, or specificity verdict.

| Case | Trace | Exit | Wall | User | System | Peak RSS |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ACTB | unset | 124 | 60.09 s | 58.28 s | 1.30 s | 932,680 KiB |
| ACTB | `1` | 124 | 60.08 s | 58.37 s | 1.31 s | 931,736 KiB |
| TP53 | unset | 124 | 60.08 s | 58.32 s | 1.33 s | 904,552 KiB |
| TP53 | `1` | 124 | 60.07 s | 58.46 s | 1.21 s | 933,100 KiB |

All four status-124 exits are GNU `timeout` deadlines after TERM, not failed
biological verdicts. The total import wall budget was 240.32 seconds, below
300 seconds. No forced-kill status or survivor was observed.

## Phase attribution

Both traced attempts emitted 50 valid
`gentle.primer_import_trace.v1` JSONL records: 27 `begin` and 23 `completed`.
Neither trace emitted `failed`, `abandoned`, or `truncated`. Parent `import`,
`reference`, and `interpretation` spans remained open because the process was
terminated at the deadline.

The material result is concordant across ACTB and TP53:

- ACTB entered `forward_hits` at 22.602560 seconds. It emitted no terminal
  record before the 60-second deadline, so at least 37.397440 seconds of the
  bounded run elapsed inside that unfinished span.
- TP53 entered `forward_hits` at 22.694039 seconds. It emitted no terminal
  record before the deadline, so at least 37.305961 seconds elapsed inside
  that unfinished span.

Before `forward_hits`, each run performed four full annotation decodes of
465,769 records. Individual decode spans took 3.565-4.145 seconds. Each run
also completed `subject_aliases` once: 5.212 seconds for ACTB and 5.260 seconds
for TP53, each with count 1,863,076. Database inspection, repeated validation,
annotation hash/read/decode, tabular parsing, and annotation lookup are all
separately observable in `trace-analysis.json`.

The forward-hit span did not complete, so no forward-hit count, reverse-hit
evaluation, amplicon join, completed import, saved state, or output summary is
available. Nested spans are not summed. These records attribute bounded wall
time; they do not provide CPU stacks or prove which function inside
`forward_hits` consumes CPU. The same-deadline baseline/traced timings do not
establish a speedup or tracing-overhead estimate.

## Finding

The previously missing phase-attribution result is now present: on both
retained real cases, the ordinary dev build reaches `forward_hits` after about
22.6 seconds and remains there for the rest of the 60-second window. Repeated
full annotation decoding and alias construction are also measurable earlier
costs. No optimization or removal of validation follows automatically from
this evidence.

