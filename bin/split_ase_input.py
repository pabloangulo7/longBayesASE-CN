#!/usr/bin/env python3
"""Validate contrasts and balance genes across independent ASE shards.

Several contrasts can be requested at once. Each one is sharded on its own,
so a shard always holds a single test and the model never sees two.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REQUIRED = {
    "sample", "ID", "H1_counts", "H2_counts", "NonHS_counts",
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--contrast", required=True,
                        help="condition_a:condition_b, or several separated by commas")
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
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
        if "group" not in source_fields:
            raise ValueError("ASE input requires a group column")
        rows = list(reader)
    fieldnames = source_fields + ["test", "groupA", "groupB"]
    available = {row["group"] for row in rows}
    manifest = []
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
            manifest.append((contrast, group_a, group_b, index, len(genes), loads[index - 1], path.name))
    with open(output / "manifest.tsv", "w") as handle:
        handle.write("contrast\tgroup_a\tgroup_b\tshard\tgenes\ttotal_counts\tfile\n")
        for row in manifest:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
