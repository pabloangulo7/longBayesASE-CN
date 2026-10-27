#!/usr/bin/env python3
"""Merge ASE shard result files with a deterministic order."""

from __future__ import annotations

import argparse
import csv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = []
    fields = None
    for path in args.inputs:
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fields is None:
                fields = reader.fieldnames
            elif reader.fieldnames != fields:
                raise ValueError(f"ASE result schema differs in {path}")
            rows.extend(reader)
    rows.sort(key=lambda row: (row.get("test", ""), row.get("ID", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["ID", "test"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
