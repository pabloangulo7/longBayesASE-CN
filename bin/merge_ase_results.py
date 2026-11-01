#!/usr/bin/env python3
"""Merge the ASE shard results and mark the fitted features with uninformative priors.

A fitted feature gets status low_prior when min(r1, r2) <= --min-prior, because
few reads of one haplotype cover a heterozygous site, so its counts barely
inform alpha and a small error in r1 or r2 becomes a large error in the bias
correction r1/r2.
"""

from __future__ import annotations

import argparse
import csv


def load_low_prior(path: str, min_prior: float) -> set[str]:
    with open(path) as handle:
        return {row["ID"] for row in csv.DictReader(handle, delimiter="\t")
                if min(float(row["H1_prior"]), float(row["H2_prior"])) <= min_prior}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--ase-input", required=True)
    parser.add_argument("--min-prior", type=float, default=0.01)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    low_prior = load_low_prior(args.ase_input, args.min_prior)
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
        if row["ID"] in low_prior and row["status"] in ("pass", "ev_bound"):
            row["status"] = "low_prior"
    rows.sort(key=lambda row: (row.get("test", ""), row.get("ID", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["ID", "test", "status"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
