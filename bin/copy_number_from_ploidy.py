#!/usr/bin/env python3
"""Expand chromosome-level sample ploidy into per-gene copy numbers."""

from __future__ import annotations

import argparse
import csv

from common import MIN_COPY_NUMBER, gene_chromosomes, parse_ploidy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.annotation_summary) as handle:
        genes = [(row["GID"], *gene_chromosomes(row)) for row in csv.DictReader(handle, delimiter="\t")]
    dna = {}
    with open(args.samplesheet) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dna.setdefault(row["dna_id"], parse_ploidy(row.get("ploidy", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "ID", "Chrom_H1", "Chrom_H2", "CN_H1", "CN_H2"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in sorted(dna):
            ploidy = dna[sample]
            for gene, chrom_h1, chrom_h2 in genes:
                # Each allele takes the ploidy of the homologue it lies on.
                h1 = max(ploidy.get(chrom_h1, (1.0, 1.0))[0], MIN_COPY_NUMBER)
                h2 = max(ploidy.get(chrom_h2, (1.0, 1.0))[1], MIN_COPY_NUMBER)
                writer.writerow({"sample": sample, "ID": gene, "Chrom_H1": chrom_h1, "Chrom_H2": chrom_h2,
                                 "CN_H1": h1, "CN_H2": h2})


if __name__ == "__main__":
    main()
