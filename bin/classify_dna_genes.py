#!/usr/bin/env python3
"""Classify DNA reads by the genes their co-optimal alignments overlap."""

from __future__ import annotations

import argparse
import csv
import itertools
import re
from collections import Counter, defaultdict


CATEGORIES = [
    "H1", "H1_multigene", "H1_multimapping", "H1_multimapping_multigene",
    "H1_multimapping_multigene_copy", "H1_multimapping_multigene_diff",
    "H2", "H2_multigene", "H2_multimapping", "H2_multimapping_multigene",
    "H2_multimapping_multigene_copy", "H2_multimapping_multigene_diff",
    "NonHS", "NonHS_multigene", "NonHS_multimapping", "NonHS_multimapping_multigene",
    "NonHS_multigene_copy", "NonHS_multimapping_multigene_copy",
    "NonHS_multimapping_multigene_diff", "NonHS_complex", "NonHS_multimapping_complex",
    "Unclassified",
]
COMPLEX_CATEGORIES = {
    "NonHS_complex", "NonHS_multimapping_complex", "H1_multimapping_multigene_diff",
    "H2_multimapping_multigene_diff", "NonHS_multimapping_multigene_diff",
}


def split_hap(gene: str) -> tuple[str, str]:
    if gene.endswith("_hap1"):
        return gene[:-5], "H1"
    if gene.endswith("_hap2"):
        return gene[:-5], "H2"
    return gene, "H1"


def clean_gene_id(gene: str) -> str:
    if gene.startswith("chr"):
        return gene
    return re.sub(r"_\d+(?=_|$)", "", gene)


def classify(info: dict, readgroup: str) -> str:
    n_genes = len(info["genes_total"])
    has_unassigned = info["has_unassigned"]
    multigene = n_genes > 1 or (has_unassigned and n_genes == 1)
    multigene_simple = len({frozenset(sub) for sub in info["genes"]}) == 1 and not has_unassigned
    multigene_copy = (
        len({frozenset(clean_gene_id(gene) for gene in sub) for sub in info["genes"]}) == 1
        and not has_unassigned
    )
    same_h1_h2 = sorted(info["genes_HP"]["H1"]) == sorted(info["genes_HP"]["H2"]) and not has_unassigned
    return (
        (
            f"{readgroup}_multigene"
            if multigene_simple
            else f"{readgroup}_multigene_copy"
            if multigene_copy
            else f"{readgroup}_complex"
            if not same_h1_h2 and "NonHS" in readgroup
            else f"{readgroup}_multigene_diff"
        )
        if multigene
        else readgroup
    ) if n_genes > 0 else "Unassigned"


def read_genes(rows) -> dict:
    """Genes of the featureCounts CORE rows of one read, one row per alignment."""
    info = {"genes": [], "genes_total": set(), "genes_HP": defaultdict(set), "has_unassigned": False}
    for row in rows:
        if len(row) != 4:
            raise ValueError(f"featureCounts CORE row must have four columns: {row}")
        _read_id, status, _n, genes_raw = row
        if status.startswith("Unassigned"):
            info["has_unassigned"] = True
            continue
        gene_list = []
        for gene_raw in genes_raw.split(","):
            gene_id, hap = split_hap(gene_raw)
            gene_list.append(gene_id)
            info["genes_total"].add(gene_id)
            info["genes_HP"][hap].add(gene_id)
        info["genes"].append(gene_list)
    return info


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--featurecounts", required=True, help="featureCounts CORE output of a BAM grouped by read")
    parser.add_argument("--readgroups", required=True, help="read groups in the same read order")
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    gene_counts = defaultdict(Counter)
    complex_reads = defaultdict(list)
    qc = Counter()
    # Both tables list the reads in the order of the BAM, so each read is
    # classified as it is read and nothing is kept per read.
    with open(args.featurecounts) as core, open(args.readgroups) as groups, \
            open(args.output_prefix + ".readgroups_genes.tsv", "w") as out:
        readgroups = csv.DictReader(groups, delimiter="\t")
        out.write("Read_ID\tGroup\n")
        for read_id, rows in itertools.groupby((row for row in csv.reader(core, delimiter="\t") if row), key=lambda row: row[0]):
            row = next(readgroups, None)
            if row is None or row["Read_ID"] != read_id:
                raise ValueError(f"read {read_id}: the featureCounts and read-group tables are not in the same read order")
            info = read_genes(rows)
            group = classify(info, row["Group"])
            out.write(f"{read_id}\t{group}\n")
            qc[group] += 1
            genes = info["genes_total"]
            if not genes:
                continue
            genes_to_count = {clean_gene_id(gene) for gene in genes} if group.endswith("_copy") else genes
            for gene in genes_to_count:
                if group in COMPLEX_CATEGORIES:
                    complex_reads[gene].append(read_id)
                gene_counts[gene][group] += 1

    with open(args.output_prefix + ".gene_counts.tsv", "w") as out:
        out.write("GID\t" + "\t".join(CATEGORIES) + "\n")
        for gene, counts in gene_counts.items():
            out.write(gene + "\t" + "\t".join(str(counts[cat]) for cat in CATEGORIES) + "\n")
    with open(args.output_prefix + ".complex_reads.tsv", "w") as out:
        out.write("GID\tReads\n")
        for gene, reads in complex_reads.items():
            out.write(f"{gene}\t{','.join(reads)}\n")
    with open(args.output_prefix + ".mapping_qc.tsv", "w") as out:
        out.write("category\treads\n")
        for category, count in sorted(qc.items()):
            out.write(f"{category}\t{count}\n")


if __name__ == "__main__":
    main()
