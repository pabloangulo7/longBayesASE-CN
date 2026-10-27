#!/usr/bin/env python3
"""Filter the ASE input and balance its features across independent shards.

Features whose mapping priors cannot inform the model are left out first:
those whose reads rarely cover a heterozygous site on one of the haplotypes
(min(r1, r2) <= --min-prior), and those whose simulated reads come out on the
wrong haplotype more often than --max-misassignment, which the model has no
term for. Each contrast is then sharded on its own, so a shard always holds a
single test and the model never sees two.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


REQUIRED = {
    "sample", "group", "ID", "H1_counts", "H2_counts", "NonHS_counts",
    "CN_H1", "CN_H2", "H1_prior", "H2_prior",
}


def parse_contrast(value: str):
    """condition_a:condition_b, exactly as the two appear in the samplesheet."""
    parts = value.split(":")
    if len(parts) != 2 or not all(part.strip() for part in parts):
        raise ValueError(f"contrast '{value}' must use condition_a:condition_b syntax")
    group_a, group_b = (part.strip() for part in parts)
    if group_a == group_b:
        raise ValueError(f"contrast '{value}' compares a condition with itself")
    return f"{group_a}_VS_{group_b}", group_a, group_b


def parse_contrasts(value: str):
    """One or more contrasts, comma separated."""
    items = [item.strip() for item in value.split(",") if item.strip()]
    if not items:
        raise ValueError("--contrast is empty")
    contrasts = [parse_contrast(item) for item in items]
    labels = [contrast[0] for contrast in contrasts]
    repeated = sorted({label for label in labels if labels.count(label) > 1})
    if repeated:
        raise ValueError(f"repeated contrast(s): {', '.join(repeated)}")
    return contrasts


def exclusion(row: dict[str, str], min_prior: float, max_misassignment: float) -> str:
    """Why a feature cannot be fitted, or an empty string if it can."""
    if min(float(row["H1_prior"]), float(row["H2_prior"])) <= min_prior:
        return f"min(r1, r2) <= {min_prior}"
    wrong = max(float(row.get("H1_misassignment") or 0), float(row.get("H2_misassignment") or 0))
    if wrong > max_misassignment:
        return f"misassignment > {max_misassignment}"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--contrast", required=True,
                        help="condition_a:condition_b, or several separated by commas")
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--min-prior", type=float, default=0.01)
    parser.add_argument("--max-misassignment", type=float, default=0.05)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--excluded", required=True, help="table of the features left out and why")
    args = parser.parse_args()
    contrasts = parse_contrasts(args.contrast)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with open(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"ASE input is missing columns: {', '.join(sorted(missing))}")
        source_fields = list(reader.fieldnames or [])
        rows = list(reader)
    excluded = {}
    for row in rows:
        reason = exclusion(row, args.min_prior, args.max_misassignment)
        if reason:
            excluded.setdefault(row["ID"], reason)
    with open(args.excluded, "w") as handle:
        handle.write("ID\treason\n")
        handle.writelines(f"{feature}\t{reason}\n" for feature, reason in sorted(excluded.items()))
    print(f"[split_ase_input] left out {len(excluded)} features; see {args.excluded}", file=sys.stderr)
    available = {row["group"] for row in rows}
    rows = [row for row in rows if row["ID"] not in excluded]
    if not rows:
        raise ValueError("no feature passed the mapping-prior filters")
    fieldnames = source_fields + ["test", "groupA", "groupB"]
    for contrast, group_a, group_b in contrasts:
        missing_groups = sorted({group_a, group_b} - available)
        if missing_groups:
            raise ValueError(
                f"contrast {contrast} references absent conditions "
                f"({', '.join(missing_groups)}); available: {', '.join(sorted(available))}")
        by_gene = defaultdict(list)
        for row in rows:
            if row["group"] not in (group_a, group_b):
                continue
            row = dict(row, test=contrast, groupA=group_a, groupB=group_b)
            by_gene[row["ID"]].append(row)
        weights = {
            gene: sum(int(float(row["H1_counts"])) + int(float(row["H2_counts"])) + int(float(row["NonHS_counts"])) for row in gene_rows)
            for gene, gene_rows in by_gene.items()
        }
        shard_count = min(max(1, args.shards), len(by_gene))
        bins = [[] for _ in range(shard_count)]
        loads = [0] * shard_count
        for gene in sorted(by_gene, key=lambda item: (-weights[item], item)):
            index = min(range(shard_count), key=lambda item: loads[item])
            bins[index].append(gene)
            loads[index] += weights[gene]
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast)
        for index, genes in enumerate(bins, 1):
            path = output / f"{safe}.shard_{index:04d}.tsv"
            with open(path, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for gene in genes:
                    writer.writerows(by_gene[gene])


if __name__ == "__main__":
    main()
