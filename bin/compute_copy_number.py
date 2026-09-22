#!/usr/bin/env python3
"""Estimate chromosome-wide haplotype copy number from ONT gDNA depth."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

from common import MIN_COPY_NUMBER, feature_parts, median, parse_ploidy


DIPLOID = 2.0


COVERAGE_RE = re.compile(
    r"^(?P<dna>.+?)[._](?P<group>HS|NonHS|HS_copy|NonHS_copy)[._]coverage\.(?:tsv|txt)$"
)


def finite_median(values) -> float:
    values = [value for value in values if math.isfinite(value)]
    return median(values) if values else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--coverage", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--qc-output", required=True)
    parser.add_argument("--min-haplotype-depth", type=float, default=10.0,
                        help="haplotype-specific depth a gene needs before its "
                             "own H1/H2 split is trusted")
    parser.add_argument("--min-gene-depth", type=float, default=1.0,
                        help="total depth a gene needs before its own copy "
                             "number is estimated at gene resolution")
    parser.add_argument("--gene-copies", choices=("all", "primary"), default="all",
                        help="sum the coverage of every copy Liftoff found, or "
                             "only the primary one")
    parser.add_argument("--cn-change-threshold", type=float, default=1.2,
                        help="fold change from the diploid baseline, in either "
                             "direction, that makes a chromosome aneuploid")
    parser.add_argument("--cn-resolution", choices=("chromosome", "gene"), default="chromosome")
    args = parser.parse_args()
    if args.cn_change_threshold <= 1:
        parser.error("--cn-change-threshold must be above 1")

    with open(args.annotation_summary) as handle:
        annotations = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row.get("GID") or not row.get("Chrom"):
                raise ValueError("annotation summary requires GID and Chrom")
            annotations[row["GID"]] = row
    if not annotations:
        raise ValueError("annotation summary contains no genes")

    expected_ploidy: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    dna_ids = set()
    with open(args.samplesheet) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dna = row["dna_id"]
            dna_ids.add(dna)
            expected_ploidy[dna].update(parse_ploidy(row.get("ploidy", "")))

    depth: dict[tuple[str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    seen_groups: dict[str, set[str]] = defaultdict(set)
    for path in args.coverage:
        match = COVERAGE_RE.match(Path(path).name)
        if not match:
            raise ValueError(f"unrecognized coverage filename: {path}")
        dna, group = match.group("dna"), match.group("group")
        if group in {"HS", "NonHS"}:
            seen_groups[dna].add(group)
        with open(path) as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not row or row[0].startswith("#"):
                    continue
                gid, hap = feature_parts(row[0])
                if gid not in annotations:
                    continue
                value = float(row[-1])
                if group == "HS" and hap:
                    depth[(dna, gid)][hap] += value
                elif group == "NonHS" and hap == "H1":
                    depth[(dna, gid)]["NonHS"] += value
                elif group == "HS_copy" and hap:
                    depth[(dna, gid)][f"{hap}_copy"] += value
                elif group == "NonHS_copy" and hap == "H1":
                    depth[(dna, gid)]["NonHS_copy"] += value

    for dna in dna_ids:
        missing = {"HS", "NonHS"} - seen_groups[dna]
        if missing:
            raise ValueError(f"DNA input {dna} is missing coverage groups: {', '.join(sorted(missing))}")

    output_rows = []
    qc_rows = []
    for dna in sorted(dna_ids):
        genes = []
        for gid, annotation in annotations.items():
            values = depth[(dna, gid)]
            h1, h2, nonhs = values["H1"], values["H2"], values["NonHS"]
            h1_copy, h2_copy, nonhs_copy = (
                values["H1_copy"], values["H2_copy"], values["NonHS_copy"]
            )
            # Liftoff reports additional copies of a gene elsewhere in the
            # assembly. Summing them measures how many copies of that sequence
            # the genome carries, which is what dosage means for expression;
            # keeping the primary one measures the canonical locus alone.
            extra = args.gene_copies == "all"
            hs = h1 + h2 + (h1_copy + h2_copy if extra else 0.0)
            total = hs + nonhs + (nonhs_copy if extra else 0.0)
            h1_depth = h1 + (h1_copy if extra else 0.0)
            genes.append(
                {
                    "ID": gid,
                    "Chrom": annotation["Chrom"],
                    "Set": annotation.get("Set", "NA"),
                    "H1_mean": h1,
                    "H2_mean": h2,
                    "NonHS_mean": nonhs,
                    "H1_copy_mean": h1_copy,
                    "H2_copy_mean": h2_copy,
                    "NonHS_copy_mean": nonhs_copy,
                    "Total_mean": total,
                    "HS_mean": hs,
                    # Below the haplotype-specific depth the split is noise, so
                    # the gene later borrows its chromosome's proportion.
                    "gene_H1_prop": h1_depth / hs if hs >= args.min_haplotype_depth else float("nan"),
                }
            )

        baseline_values = [
            row["Total_mean"]
            for row in genes
            if row["Set"] == "Singleton"
            and (
                row["Chrom"] not in expected_ploidy[dna]
                or sum(expected_ploidy[dna][row["Chrom"]]) == DIPLOID
            )
        ]
        baseline = finite_median(baseline_values)
        if not math.isfinite(baseline) or baseline <= 0:
            raise ValueError(f"cannot estimate a positive diploid coverage baseline for {dna}")

        singletons: dict[str, list[dict]] = defaultdict(list)
        for row in genes:
            if row["Set"] == "Singleton":
                singletons[row["Chrom"]].append(row)
        global_h1_prop = finite_median(row["gene_H1_prop"] for row in genes if row["Set"] == "Singleton")
        if not math.isfinite(global_h1_prop):
            global_h1_prop = 0.5

        chrom_stats = {}
        for chrom in sorted({row["Chrom"] for row in genes}):
            chrom_genes = singletons.get(chrom, [])
            # A declared ploidy only keeps the chromosome out of the baseline,
            # which happens above. The depth itself is always what decides the
            # copy number here, because it measures what the cells actually
            # carry: a subclonal or partial gain is not the declared integer.
            chrom_cn = finite_median(
                DIPLOID * row["Total_mean"] / baseline for row in chrom_genes
            )
            if not math.isfinite(chrom_cn):
                chrom_cn = DIPLOID
            ratio = chrom_cn / DIPLOID
            # Symmetric: a loss is as much a departure from diploidy as a gain.
            changed = ratio >= args.cn_change_threshold or ratio <= 1.0 / args.cn_change_threshold
            h1_prop = finite_median(row["gene_H1_prop"] for row in chrom_genes)
            if not math.isfinite(h1_prop):
                h1_prop = global_h1_prop
            h1_prop = min(max(h1_prop, 0.0), 1.0)
            if changed:
                cn_h1 = max(chrom_cn * h1_prop, MIN_COPY_NUMBER)
                cn_h2 = max(chrom_cn * (1.0 - h1_prop), MIN_COPY_NUMBER)
            else:
                cn_h1 = cn_h2 = 1.0
            chrom_stats[chrom] = (chrom_cn, ratio, changed, h1_prop, cn_h1, cn_h2)
            qc_rows.append(
                {
                    "sample": dna,
                    "Chrom": chrom,
                    "baseline_depth": baseline,
                    "baseline_genes": len(baseline_values),
                    "chromosome_singleton_genes": len(chrom_genes),
                    "haplotype_informative_genes": sum(
                        math.isfinite(row["gene_H1_prop"]) for row in chrom_genes
                    ),
                    "expected_ploidy": chrom in expected_ploidy[dna],
                    "Chrom_CN": chrom_cn,
                    "ratio": ratio,
                    "changed": changed,
                    "direction": "gain" if ratio >= args.cn_change_threshold else
                                 "loss" if ratio <= 1.0 / args.cn_change_threshold else "none",
                    "H1_prop": h1_prop,
                    "H2_prop": 1.0 - h1_prop,
                    "CN_H1": cn_h1,
                    "CN_H2": cn_h2,
                }
            )

        for row in genes:
            chrom_cn, ratio, changed, h1_prop, cn_h1, cn_h2 = chrom_stats[row["Chrom"]]
            if args.cn_resolution == "gene" and row["Total_mean"] >= args.min_gene_depth:
                # A gene with no coverage at all cannot tell a deletion from a
                # region its reads never reach, so it keeps its chromosome's
                # value instead of being called homozygously deleted.
                gene_prop = row["gene_H1_prop"]
                if not math.isfinite(gene_prop):
                    gene_prop = h1_prop
                cn_total = DIPLOID * row["Total_mean"] / baseline
                cn_h1 = max(cn_total * gene_prop, MIN_COPY_NUMBER)
                cn_h2 = max(cn_total * (1.0 - gene_prop), MIN_COPY_NUMBER)
            output_rows.append(
                {
                    "sample": dna,
                    "ID": row["ID"],
                    "Chrom": row["Chrom"],
                    "CN_H1": cn_h1,
                    "CN_H2": cn_h2,
                }
            )

    for path, rows in ((args.output, output_rows), (args.qc_output, qc_rows)):
        if not rows:
            raise ValueError(f"no rows available for {path}")
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
