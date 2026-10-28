#!/usr/bin/env python3
"""Merge the ASE shard results and flag the features with uninformative priors.

Every feature is fitted; filter_flag is low_prior when min(r1, r2) <= --min-prior,
because few reads of one haplotype cover a heterozygous site, so its counts
barely inform alpha and a small error in r1 or r2 becomes a large error in the
bias correction r1/r2. pass otherwise.
"""

from __future__ import annotations

import argparse
import csv


def load_flags(path: str, min_prior: float) -> dict[str, str]:
    with open(path) as handle:
        return {
            row["ID"]: "low_prior" if min(float(row["H1_prior"]), float(row["H2_prior"])) <= min_prior else "pass"
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--ase-input", required=True)
    parser.add_argument("--min-prior", type=float, default=0.01)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    flags = load_flags(args.ase_input, args.min_prior)
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
    for row in rows:
        row["filter_flag"] = flags[row["ID"]]
    rows.sort(key=lambda row: (row.get("test", ""), row.get("ID", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*(fields or ["ID", "test"]), "filter_flag"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
