#!/usr/bin/env python3
"""Classify RNA reads as H1, H2 or NonHS and count them per gene or isoform.

A read is locked to the first assignment whose probability reaches
--resolve-threshold; otherwise it keeps every feature it was assigned to.

--min-probability drops the assignments below it before a read is classified,
so a read whose evidence is concentrated on one feature is no longer discarded
for the trace probabilities left on the others. A read with nothing above the
threshold is counted as Unclassified and attributed to no feature.

The assignment table lists all the rows of a read together (it is written from
Oarfish one read at a time), so reads are classified as they stream past and
memory does not grow with library depth. A read whose rows are split into two
blocks stops the run instead of being counted twice.
"""

from __future__ import annotations

import argparse
import hashlib
from array import array
from collections import Counter, defaultdict

import numpy as np


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

H1, H2 = 0, 1


def split_hap(feature: str) -> tuple[str, int]:
    if feature.endswith("_hap1"):
        return feature[:-5], H1
    if feature.endswith("_hap2"):
        return feature[:-5], H2
    raise ValueError(f"feature lacks _hap1/_hap2 suffix: {feature}")


def read_key(read_id: str) -> int:
    return int.from_bytes(hashlib.blake2b(read_id.encode(), digest_size=8).digest(), "little", signed=True)


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
    parser.add_argument("--min-probability", type=float, default=0.0)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    if not 0 < args.resolve_threshold <= 1:
        raise ValueError("--resolve-threshold must be in (0,1]")
    if not 0 <= args.min_probability < args.resolve_threshold:
        raise ValueError("--min-probability must be in [0,--resolve-threshold)")

    gene_counts: dict[str, Counter] = defaultdict(Counter)
    complex_reads: dict[str, list[str]] = defaultdict(list)
    qc: Counter = Counter()
    features: dict[str, str] = {}
    use_gene = args.mode == "genes"

    with open(args.assignments) as handle, open(args.output_prefix + ".readgroups.tsv", "w") as groups_out:
        groups_out.write("Read_ID\tGroup\n")

        def finish(read_id: str, h1: list[str], h2: list[str]) -> None:
            genes = set(h1)
            genes.update(h2)
            group = classify(h1, h2, genes) if genes else "Unclassified"
            groups_out.write(f"{read_id}\t{group}\n")
            qc[group] += 1
            complex_group = group in COMPLEX_CATEGORIES
            for gene in genes:
                if complex_group:
                    complex_reads[gene].append(read_id)
                gene_counts[gene][group] += 1

        current = None
        record: list = [[], []]
        locked = False
        # 8 bytes per read instead of the read names, to check at the end that
        # no read was split into two blocks.
        seen = array("q")
        for line_number, line in enumerate(handle, 1):
            row = line.rstrip("\n").split("\t")
            if row == [""] or row[0] == "read_id":
                continue
            if len(row) != 4:
                raise ValueError(f"line {line_number}: assignment row must have four columns: {row}")
            read_id, tx_id, gene_id, probability = row
            if read_id != current:
                if current is not None:
                    finish(current, record[H1], record[H2])
                # Every read is registered, including one left with no
                # assignment above --min-probability, so the check below still
                # sees a read whose rows come back in a later block.
                seen.append(read_key(read_id))
                current, record, locked = read_id, [[], []], False
            elif locked:
                continue
            value = float(probability)
            if value < args.min_probability:
                continue
            feature, hap = split_hap(gene_id if use_gene else tx_id)
            # One shared string per feature instead of one copy per read.
            feature = features.setdefault(feature, feature)
            if value >= args.resolve_threshold:
                record = [[], []]
                locked = True
            record[hap].append(feature)
        if current is not None:
            finish(current, record[H1], record[H2])

    keys = np.frombuffer(seen, dtype=np.int64)
    keys.sort()
    if keys.size > 1 and np.any(keys[1:] == keys[:-1]):
        raise ValueError(f"{args.assignments}: the rows of at least one read are not contiguous")

    with open(args.output_prefix + ".counts.tsv", "w") as out:
        out.write("ID\t" + "\t".join(CATEGORIES) + "\n")
        # Sorted so that two runs over the same input give the same file.
        for gene, counts in sorted(gene_counts.items()):
            out.write(gene + "\t" + "\t".join(str(counts[cat]) for cat in CATEGORIES) + "\n")
    with open(args.output_prefix + ".complex_reads.tsv", "w") as out:
        out.write("ID\tReads\n")
        for gene, reads in sorted(complex_reads.items()):
            out.write(f"{gene}\t{','.join(reads)}\n")
    with open(args.output_prefix + ".qc.tsv", "w") as out:
        out.write("metric\treads\n")
        for metric, value in sorted(qc.items()):
            out.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
