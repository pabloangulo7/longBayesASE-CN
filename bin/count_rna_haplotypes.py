#!/usr/bin/env python3
"""Classify RNA reads as H1, H2 or NonHS and count them per gene or isoform.

A read is locked to the first assignment whose probability reaches
--resolve-threshold; otherwise it keeps every feature it was assigned to.

A deep cDNA library has tens of millions of reads and every one of them is held
in memory until the file has been read, so the per-read record is kept as small
as it can be: two lists of features, one per haplotype, and the lock flag. The
set of features a read touches is the union of those two lists, so it is derived
rather than stored, and the feature strings are shared instead of copied.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict


CATEGORIES = [
    "H1", "H1_multimapping", "H1_multimapping_multigene",
    "H2", "H2_multimapping", "H2_multimapping_multigene",
    "NonHS", "NonHS_multimapping", "NonHS_multimapping_multigene",
    "NonHS_complex", "NonHS_multimapping_complex", "Unclassified",
]
COMPLEX_CATEGORIES = {
    "H1_multimapping_multigene", "H2_multimapping_multigene",
    "NonHS_multimapping_multigene", "NonHS_complex", "NonHS_multimapping_complex",
}

H1, H2, LOCKED = 0, 1, 2


def split_hap(feature: str) -> tuple[str, int]:
    if feature.endswith("_hap1"):
        return feature[:-5], H1
    if feature.endswith("_hap2"):
        return feature[:-5], H2
    raise ValueError(f"feature lacks _hap1/_hap2 suffix: {feature}")


def classify(h1: list[str], h2: list[str], genes: set[str]) -> str:
    haplotype = "NonHS" if h1 and h2 else ("H1" if h1 else "H2")
    multimapping = len(h1) > 1 or len(h2) > 1
    multigene = len(genes) > 1
    same_h1_h2 = sorted(set(h1)) == sorted(set(h2))
    return (
        ("NonHS_multimapping_complex" if multimapping else "NonHS_complex")
        if haplotype == "NonHS" and not same_h1_h2
        else (f"{haplotype}_multimapping_multigene" if multigene else f"{haplotype}_multimapping")
        if multimapping
        else haplotype
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--mode", choices=["genes", "isoforms"], required=True)
    parser.add_argument("--resolve-threshold", type=float, default=0.90)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    if not 0 < args.resolve_threshold <= 1:
        raise ValueError("--resolve-threshold must be in (0,1]")

    read_info: dict[str, list] = {}
    features: dict[str, str] = {}
    with open(args.assignments) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0] == "read_id":
                continue
            if len(row) != 4:
                raise ValueError(f"assignment row must have four columns: {row}")
            read_id, tx_id, gene_id, probability = row
            record = read_info.get(read_id)
            if record is None:
                record = read_info[read_id] = [[], [], False]
            elif record[LOCKED]:
                continue
            feature, hap = split_hap(gene_id if args.mode == "genes" else tx_id)
            # One shared string per feature instead of one copy per read.
            feature = features.setdefault(feature, feature)
            if float(probability) >= args.resolve_threshold:
                record[H1], record[H2], record[LOCKED] = [], [], True
                record[hap].append(feature)
                continue
            record[hap].append(feature)

    gene_counts = defaultdict(Counter)
    complex_reads = defaultdict(list)
    qc = Counter()
    # Classifying, writing the read groups and counting in one pass keeps the
    # reads out of a second structure.
    with open(args.output_prefix + ".readgroups.tsv", "w") as out:
        out.write("Read_ID\tGroup\n")
        for read_id, (h1, h2, _locked) in read_info.items():
            genes = set(h1)
            genes.update(h2)
            group = classify(h1, h2, genes)
            out.write(f"{read_id}\t{group}\n")
            qc[group] += 1
            complex_group = group in COMPLEX_CATEGORIES
            for gene in genes:
                if complex_group:
                    complex_reads[gene].append(read_id)
                gene_counts[gene][group] += 1

    with open(args.output_prefix + ".counts.tsv", "w") as out:
        out.write("ID\t" + "\t".join(CATEGORIES) + "\n")
        for gene, counts in gene_counts.items():
            out.write(gene + "\t" + "\t".join(str(counts[cat]) for cat in CATEGORIES) + "\n")
    with open(args.output_prefix + ".complex_reads.tsv", "w") as out:
        out.write("ID\tReads\n")
        for gene, reads in complex_reads.items():
            out.write(f"{gene}\t{','.join(reads)}\n")
    with open(args.output_prefix + ".qc.tsv", "w") as out:
        out.write("metric\treads\n")
        for metric, value in sorted(qc.items()):
            out.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
