#!/usr/bin/env python3
"""Merge the ASE shard results and add the Bayesian false sign rates.

For every gene, diffAI_lfsr = diffAI_pvalue / 2 is the posterior probability
that the sign of the difference between the groups is wrong (local false sign
rate). Within each test, diffAI_bfdr of a gene is the mean lfsr of all the
genes at least as significant: the expected proportion of wrong directions in
the list that stops at that gene. Both are reported, not used to filter.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict

NEW_COLUMNS = ["diffAI_lfsr", "diffAI_bfdr"]


def as_float(value: str | None) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def add_bayesian_fdr(rows: list[dict[str, str]]) -> None:
    by_test: dict[str, list[tuple[float, dict[str, str]]]] = defaultdict(list)
    for row in rows:
        lfsr = as_float(row.get("diffAI_pvalue")) / 2
        if math.isfinite(lfsr):
            by_test[row["test"]].append((lfsr, row))
        else:
            row.update(dict.fromkeys(NEW_COLUMNS, "NA"))
    for genes in by_test.values():
        genes.sort(key=lambda item: item[0])
        cumulative = 0.0
        bfdr = []
        for rank, (lfsr, _row) in enumerate(genes, start=1):
            cumulative += lfsr
            bfdr.append(cumulative / rank)
        # Genes tied on lfsr share the rate of the whole tied group.
        for index in range(len(genes) - 2, -1, -1):
            if genes[index][0] == genes[index + 1][0]:
                bfdr[index] = bfdr[index + 1]
        for (lfsr, row), rate in zip(genes, bfdr):
            row["diffAI_lfsr"] = f"{lfsr:.6g}"
            row["diffAI_bfdr"] = f"{rate:.6g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
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
    fields = list(fields or ["ID", "test"])
    if "diffAI_pvalue" in fields:
        add_bayesian_fdr(rows)
        position = fields.index("diffAI_pvalue") + 1
        fields[position:position] = NEW_COLUMNS
    rows.sort(key=lambda row: (row.get("test", ""), row.get("ID", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
