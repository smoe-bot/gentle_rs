# TSS Tutorial Lifecycle Acceptance — 2026-09-19

## Scope And Verdicts

This report covers the public synthetic `tss_collection_gui` tutorial only. It
does not cover TP73/DeltaNp73, private data, the inner agent, active-promoter
biology or broader release readiness.

| Scope | Verdict | Boundary |
| --- | --- | --- |
| Deterministic source/tests | Pass | Tutorial, runner and focused Rust contracts only |
| Automated Linux GUI | Pass | Two fresh, network-isolated 18-step X11 runs |
| Manual lifecycle | Pass | Synthetic save/restart, visual orientation, deliberate stale-member edit, Forget and Undo |
| Biological/real-data acceptance | Not performed | TP73/DeltaNp73 remains a separate private-data exercise |

Baseline: `c2dead07aa021fb5b3f84274ba088d7bb923a05f`.

Runtime candidate after reproduced fixes:
`c3c06da35acd01814cd0630d278d8a1660888d1a`. The PR may contain a later
report-only commit; runtime claims remain bound to the exact candidate and
binary hashes below.

## Frozen Environment

- Debian forky/sid, Linux `7.1.3+deb14-amd64`, x86-64
- `rustc 1.95.0`, `cargo 1.95.0`, Python `3.14.7`, Git `2.55.0`
- Xvfb at `1600x1000x24`, Openbox EWMH window manager, `xdotool`, `xdpyinfo`,
  `xprop`, `xwininfo` and `scrot 2.0.0`
- fresh `HOME`, XDG config/cache/data and temporary directories per automated run
- network isolation enforced with a new Linux user/network namespace; each
  ledger records parent and child namespace identities
- `Cargo.lock` SHA-256:
  `67a96b07e6e2d5beb5b0fd1eb2ba57e1c11b6d1a4a355dab5bce6ca830a25236`

Candidate binaries were built together with
`cargo build --locked --features gui-test-support --bin gentle --bin gentle_cli --bin gentle_examples_docs`:

| Binary | SHA-256 |
| --- | --- |
| `gentle` | `5736f693571afd9e87fba3d4bf2bad2dd98b62088d747eeef53c86b869321ba1` |
| `gentle_cli` | `8c606f205001e87c3bf7f95cbf4caffd80303fba14e8537e0405e4e8b98465c8` |
| `gentle_examples_docs` | `7252fc77f9ec7455351625b90de55a28ce2f324dc8d1110a3dca82e2ddd34644` |

## Reproduced Defects And Narrow Fixes

1. The native DNA viewer inherited an `800x600` viewport and a content-sized
   toolbar panel. Wrapped `TFBS scan` controls were registered but outside the
   visible scroll area. The native viewer now starts at `1200x860`, has an
   `820x520` minimum and a bounded, user-resizable toolbar height.
2. Preview/materialization cloned the real engine snapshot on Rust's default
   2 MiB worker stack and aborted, although Cargo tests inherited the repository
   16 MiB setting. The named TSS worker now uses a bounded 16 MiB stack and
   reports worker-creation failure instead of aborting.
3. A collection-open request originated in a child viewport, while its poller
   ran only in the idle root viewport. The child now requests a root repaint
   after queueing; validation and window reuse remain application-owned.

Each change has a focused Rust regression. No sequence, TSS grouping,
approval, stale-member or scientific-verification rule was weakened.

## Automated GUI Checkpoints

Both runs started from a generated starter project, never the completed oracle.
Every row passed twice. Run A ledger SHA-256 is
`077c36161761845206552a7bf60d979bf54719b22997c42d3a3326b39da1704a`;
run B is `5f513e8e5a02651c8366415e93d540df8874ae4f33b745018c7d83e032884994`.

| Checkpoint | Expected and observed result | Verdict | Evidence |
| --- | --- | --- | --- |
| `open_locus` | Source locus opened by ordinary X11 input | Pass | A/B ledger step 1 |
| `open_tfbs_menu` | Wrapped TFBS control remained visible and enabled | Pass | A/B step 2 screenshot/snapshot |
| `open_tss_workspace` | TSS workspace opened for the source subject | Pass | A/B step 3 |
| `set_gene` | Gene field became `TOY` without mutation | Pass | A/B step 4 |
| `preview` | Three exact available starts reported | Pass | A/B step 5 typed inventory |
| `select_starts` | All available starts selected | Pass | A/B step 6 |
| `approve_derivation` | Three approved windows materialized | Pass | A/B step 7 state/oracle comparison |
| `refresh_registry` | `tss_windows` was readable but not silently validated | Pass | A/B step 8 |
| `validate_collection` | Three 701-bp members; shared membership, strands and local TSS 501 matched oracle | Pass | A/B step 9 typed report |
| `open_members` | Exactly source plus three member viewers | Pass | A/B step 10 exact subject scopes |
| `reopen_members` | Reopening reused the same four viewers | Pass | A/B step 11 exact subject scopes |
| `request_forget` | Named confirmation opened without mutation | Pass | A/B step 12 before/after facts |
| `cancel_forget` | Cancel retained collection, sequences and viewers | Pass | A/B step 13 |
| `request_forget_again` | Confirmation could be requested again | Pass | A/B step 14 |
| `confirm_forget` | Only registry metadata disappeared | Pass | A/B step 15 state and window facts |
| `open_edit_menu` | Root Edit menu exposed Undo | Pass | A/B step 16 |
| `undo_forget` | Undo restored the collection registry | Pass | A/B step 17 state facts |
| `inspect_restored_collection` | Fresh validation again matched the oracle | Pass | A/B step 18 typed report |

Screenshots establish visible interaction only. Scientific pass rows above are
backed independently by typed reports, persisted state and exact oracle sequence
content.

## Manual Lifecycle

The manual run used a disposable copy of run A's intact saved project.

| Check | Observed result | Verdict | Evidence |
| --- | --- | --- | --- |
| Rendered geometry | Minus member displayed 701 bp, genomic `toy_chr:1300..2000`, strand `-`, transcript orientation and local TSS 501 | Pass | `minus-viewer.png` SHA-256 `9f01a8050bf07fa1b3c6aa98bd7960511d17fa391c37616ec3576e2f27605b1a` |
| Save/close/restart | Reopened report retained three members and membership fingerprint `sha256:f8f6ee82…773368f7` | Pass | `reopened-collection.json` SHA-256 `6af5745cfc9d31e86316a3ae6c047d0f33bd364a563e985b9ba51306fbe6129b` |
| Deliberate damage | Feature Editor previewed and applied `misc_feature 501..501 -> 502..502` | Pass | `edit-preview.png` SHA-256 `bd46a01342a8f68668ba8b22dd493af296e496a51af6cd370930430a11c3d1e9` |
| Stale inspection | `promoters tss-collection tss_windows` returned `InvalidInput` | Pass | stderr SHA-256 `fb393804d2c20823f480abff96a617060a9b004829e0b31c4251b6096dbcd8c5` |
| Stale collection scan | TFBS collection scan returned the same stale-member `InvalidInput`, not an empty report | Pass | retained command receipt |
| Forget | Registry list became empty; all four sequences remained | Pass | state SHA-256 `97e4e600…f97ab441`, list SHA-256 `bb7673f1…59e6685` |
| Undo after Forget | Registry returned as not checked; marker remained local 502; fresh inspection still rejected it | Pass | list SHA-256 `d4ffc5d6…aad403`, rejection SHA-256 `fb393804…cd8c5` |

## Deterministic Checks

Before live acceptance, the baseline passed 32 Python runner/publication/
checkout tests, 15 focused Rust tests, catalog 58/58, manifest/tutorial 29/29,
locked Cargo check, formatting and whitespace checks. Final-revision reruns are
recorded in the PR handoff and external evidence manifest.

The direct Rust test binary requires the repository's configured
`RUST_MIN_STACK=16777216`; running that binary outside Cargo without the setting
is not an equivalent test environment.

## Evidence And Remaining Boundaries

Raw PNGs, semantic snapshots, typed receipts, project copies, command output and
logs are retained outside Git in a checksum-manifest bundle. Only this small,
publication-safe synthetic summary is committed. No private input is present.

Two pre-existing tutorial human-review staleness warnings remain unrelated to
this chapter. No real-data or inner-agent check was run, no production report
was regenerated, and no release was tagged.
