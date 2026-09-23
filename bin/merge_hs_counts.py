#!/usr/bin/env python3
"""Merge per-library haplotype-specific counts into one table.

The same code serves the gene and the transcript tables; only the input file
names differ, and both come out with the columns sample, ID, H1, H2, NonHS.
Each column adds the reads confined to one feature (for example H1 and
H1_multimapping); multigene and complex reads stay out of the model.
"""

from __future__ import annotations

import argparse
import csv

from common import files_by_sample


def total(row: dict[str, str], name: str) -> int:
    return int(float(row.get(name) or 0)) + int(float(row.get(f"{name}_multimapping") or 0))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    by_sample = files_by_sample(args.counts, r"\.(?:gene|transcript)_HS_counts\.tsv")

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
