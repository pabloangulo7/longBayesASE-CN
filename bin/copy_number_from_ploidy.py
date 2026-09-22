#!/usr/bin/env python3
"""Expand chromosome-level sample ploidy into per-gene copy numbers."""

from __future__ import annotations

import argparse
import csv

from common import MIN_COPY_NUMBER, parse_ploidy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.annotation_summary) as handle:
        genes = [(row["GID"], row["Chrom"]) for row in csv.DictReader(handle, delimiter="\t")]
    dna = {}
    with open(args.samplesheet) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dna.setdefault(row["dna_id"], parse_ploidy(row.get("ploidy", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "ID", "Chrom", "CN_H1", "CN_H2"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in sorted(dna):
            for gene, chromosome in genes:
                h1, h2 = (max(value, MIN_COPY_NUMBER) for value in dna[sample].get(chromosome, (1.0, 1.0)))
                writer.writerow({"sample": sample, "ID": gene, "Chrom": chromosome, "CN_H1": h1, "CN_H2": h2})


if __name__ == "__main__":
    main()
