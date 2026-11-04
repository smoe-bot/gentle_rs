#!/usr/bin/env python3
"""Summarize bounded GENtle primer-import JSONL traces without summing nested spans."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def analyze(path: Path) -> dict:
    records = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        record = json.loads(line)
        if record.get("schema") != "gentle.primer_import_trace.v1":
            raise ValueError(f"unexpected schema at {path}:{number}")
        records.append(record)

    begins = {(r["import"], r["span"]): r for r in records if r["kind"] == "begin"}
    terminals = {
        (r["import"], r["span"]): r
        for r in records
        if r["kind"] in {"completed", "failed", "abandoned"}
    }
    completed = [r for r in records if r["kind"] == "completed"]
    unmatched = [begins[key] for key in begins.keys() - terminals.keys()]
    last_elapsed = max((r["elapsed_us"] for r in records), default=0)

    return {
        "path": str(path),
        "record_count": len(records),
        "schema_counts": dict(Counter(r["schema"] for r in records)),
        "kind_counts": dict(Counter(r["kind"] for r in records)),
        "phase_counts": dict(Counter(r["phase"] for r in records)),
        "completed_spans": [
            {
                "import": r["import"],
                "span": r["span"],
                "phase": r["phase"],
                "duration_us": r["duration_us"],
                **({"count": r["count"]} if "count" in r else {}),
            }
            for r in completed
        ],
        "unmatched_begins": [
            {
                "import": r["import"],
                "span": r["span"],
                "phase": r["phase"],
                "elapsed_us": r["elapsed_us"],
                "minimum_unfinished_us_at_deadline_60s": max(
                    0, 60_000_000 - r["elapsed_us"]
                ),
            }
            for r in sorted(unmatched, key=lambda item: item["span"])
        ],
        "truncated": any(r["kind"] == "truncated" for r in records),
        "last_record_elapsed_us": last_elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = {"traces": [analyze(path) for path in args.paths]}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
