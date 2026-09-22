#!/usr/bin/env python3
"""Merge per-library RNA counts into one unified table.

The same code serves the gene and the isoform tables; only the input file
names differ, and both come out with the columns sample, ID, H1, H2, NonHS.
"""

from __future__ import annotations

import argparse
import csv
import re


def total(row: dict[str, str], name: str) -> int:
    return int(float(row.get(name) or 0)) + int(float(row.get(f"{name}_multimapping") or 0))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pattern = re.compile(r"^(?P<sample>.+)\.(?:genes|isoforms)\.counts\.tsv$")
    by_sample: dict[str, str] = {}
    for path in args.counts:
        name = path.rsplit("/", 1)[-1]
        match = pattern.match(name)
        if not match:
            raise ValueError(f"cannot identify sample from RNA count filename: {name}")
        sample = match.group("sample")
        if sample in by_sample:
            raise ValueError(f"two count files for sample {sample}: {by_sample[sample]} and {path}")
        by_sample[sample] = path

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "ID", "H1", "H2", "NonHS"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Sorted by sample, then by feature, one library in memory at a time.
        for sample in sorted(by_sample):
            with open(by_sample[sample]) as source:
                rows = [
                    {
                        "sample": sample, "ID": row["ID"],
                        "H1": total(row, "H1"), "H2": total(row, "H2"),
                        "NonHS": total(row, "NonHS"),
                    }
                    for row in csv.DictReader(source, delimiter="\t")
                ]
            rows.sort(key=lambda row: row["ID"])
            writer.writerows(rows)


if __name__ == "__main__":
    main()
