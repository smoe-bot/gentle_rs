# Tutorial Coverage and Screenshot Backlog

Updated 2026-11-04 at source `6dd18fb428790c1d8817b27bc330ede04c6fc67e`.
This is the current working list; `coverage_audit.md` remains a historical
snapshot. Catalog presence, executable replay, native interaction, screenshot
quality and human scientific review are separate states.

The 2026-11-04 catalog review counted 70 entries: 30 executable manifest
chapters, 6 typed GUI contracts and 40 hand-written/reference pages (these
categories overlap). Thirty-six entries had no catalogued graphic and 20 were
still marked unreviewed. Structural regularity is strongest in the generated
chapters; a green replay does not by itself make their prose easy to follow.

## Screenshot queue

Use only synthetic or public inputs. Retain raw captures and provenance
sidecars; label crops and annotations as derived teaching images.

1. **10.03 measured gel sizing:** native lane navigation, calibrated and null
   out-of-calibration results, and GUI save/reopen.
2. **FLNA draft:** normal-size native overview/detail and orientation. The
   2400-px render and human 24-bp teaching-exon choice are already retained;
   historical annotation-release provenance remains missing.
3. **04.09 primer specificity:** native Shell/report presentation. Defer the
   full real-reference import until the measured post-BLAST performance issue
   has a verified fix.
4. **06.07 synthetic sequence design:** preview, exact proposed output,
   created DNA and same-session undo/redo.
5. **03.05 and 10.02 handoff points:** review-before-apply and sequence-
   confirmation save/reopen, for the new Synthetic Biology learning route.

Completed on 2026-11-04: **06.08 CD44 genome/transcript/protein comparison**
now retains exact-binary captures of command-palette discovery, the canonical
742-aa product and the 139-aa CD44RC product, all from one saved public locus.

Completed on 2026-11-04: **08.16 imported TSS evidence** now retains separate
complete-window and partial-query lanes plus an imported-hit detail with
strand-aware coordinates, signed score units and source hashes. Windows/macOS
and automated GUI acceptance remain pending.

Completed on 2026-11-04: **07.04 TP73 microarray evidence** now retains a
checksum-bound saved-project view of both synthetic-statistics contrast groups,
one selected vendor-derived feature with provenance and the native nearest-
feature navigation menu. Cross-platform GUI acceptance, genuine CEL replay and
independent domain interpretation remain pending.

## Human-clarity queue

Review one page at a time and prefer orientation, a visible result and a
recovery hint over additional background prose.

1. **10.04 targeted Nanopore panel:** split or layer the roughly 10,000-word
   page; provide a short successful path before the audit/reference material.
2. **10.03 measured gel sizing:** move the long appendices behind the primary
   GUI/CLI walkthrough and add visual checkpoints.
3. **09.03 protein-expression handoff:** explain the retained evidence status
   before presenting the four screenshots; distinguish teaching capture from
   accepted runtime parity.
4. **08.04 and 08.02 VKORC1:** add screenshots at the first native decision
   points. Replace machine-specific absolute paths with repository-relative
   links.
5. **07.04 microarray evidence:** native captures are now present; ask a domain
   reader to check the evidence interpretation against a genuine CEL replay.

Do not reopen the completed 08.18 repeat-cohort run, tutorial-CI registration
or release-link repair merely to obtain more images.

## Under-taught product areas

Address these as separate reviewable commits rather than one large tutorial
rewrite.

1. **Synthetic Biology integration — addressed first.** The new
   [Design, Build, Track and Verify route](./03-09_synthetic_biology_design_build_verify.md)
   connects existing unit tutorials while making their fixture and evidence
   boundaries explicit.
2. **Assembly-method selection.** Compare Gibson with shipped Golden Gate,
   Gateway, TOPO, TA/GC, In-Fusion and NEBuilder planning routes. Start with a
   bounded synthetic Golden Gate example; do not imply every preset has equal
   native or experimental acceptance.
3. **RNA structure.** The ssRNA structure interface needs a beginner exercise
   with a small public/synthetic sequence, parameter visibility and a clear
   distinction between prediction and measured structure.
4. **microRNA target scanning.** Add a compact walkthrough for target scan,
   coordinate interpretation and uncertainty; avoid expression or repression
   claims without external evidence.
5. **Guide design end to end.** Existing generated guide-filter and oligo-
   export chapters need a guided route from target/constraints through review
   and export. Do not present the current checks as a complete off-target or
   efficacy assessment.
6. **General restriction/ligation and physical sample tracking.** Extend the
   minimal digest/ligation replay and Gibson-specific arrangement/rack lessons
   into a method-neutral beginner route only after the higher-priority gaps.

Strong existing coverage in sequence inspection, Gibson planning, PCR/qPCR,
genome/transcript context, regulatory evidence, external handoffs and Sanger
confirmation should be linked rather than duplicated.
