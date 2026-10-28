#!/usr/bin/env python3
"""Estimate the mapping probabilities r1 and r2 of every feature from simulations.

r1 is the fraction of the reads simulated from H1 that are classified H1, and r2
the same for H2. Features with fewer than --min-reads simulated reads on a
haplotype get no prior.

Besides the priors, two QC tables describe the simulation without affecting the
analysis: per feature, the reads simulated from each haplotype and the fraction
classified on the wrong one (misassignment), and the summary over all features.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys


def read_counts(path: str) -> dict[str, tuple[float, float, float]]:
    """H1, H2 and NonHS reads of every feature, multimapping reads included."""
    counts = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row["ID"]] = tuple(
                float(row.get(name) or 0) + float(row.get(f"{name}_multimapping") or 0)
                for name in ("H1", "H2", "NonHS")
            )
    return counts


def write_tsv(path: str, header: list[str], rows: list[list[object]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hap1-counts", required=True)
    parser.add_argument("--hap2-counts", required=True)
    parser.add_argument("--min-reads", type=int, default=1,
                        help="simulated reads each haplotype needs for the feature to get a prior")
    parser.add_argument("--epsilon", type=float, default=0.001)
    parser.add_argument("--output", required=True, help="mapping priors: ID, H1_prior, H2_prior")
    parser.add_argument("--misassignment", required=True, help="per-feature QC of the simulation")
    parser.add_argument("--stats", required=True, help="summary of the simulation")
    args = parser.parse_args()
    from_h1, from_h2 = read_counts(args.hap1_counts), read_counts(args.hap2_counts)
    features = sorted(set(from_h1) & set(from_h2))

    def clip(value: float) -> float:
        return min(max(value, args.epsilon), 1 - args.epsilon)

    priors, qc = [], []
    totals = {"H1": [0.0, 0.0, 0.0], "H2": [0.0, 0.0, 0.0]}
    for feature in features:
        h1 = from_h1[feature]   # reads simulated from H1 called H1, H2, NonHS
        h2 = from_h2[feature]   # reads simulated from H2 called H1, H2, NonHS
        n1, n2 = sum(h1), sum(h2)
        if min(n1, n2) < max(args.min_reads, 1):
            continue
        for haplotype, calls in (("H1", h1), ("H2", h2)):
            totals[haplotype] = [total + value for total, value in zip(totals[haplotype], calls)]
        priors.append([feature, f"{clip(h1[0] / n1):.8g}", f"{clip(h2[1] / n2):.8g}"])
        qc.append([feature, f"{n1:.0f}", f"{n2:.0f}", f"{h1[1] / n1:.8g}", f"{h2[0] / n2:.8g}"])

    write_tsv(args.output, ["ID", "H1_prior", "H2_prior"], priors)
    write_tsv(args.misassignment, ["ID", "n_hap1", "n_hap2", "H1_misassignment", "H2_misassignment"], qc)

    def median(rows: list[list[object]], column: int) -> str:
        return f"{statistics.median(float(row[column]) for row in rows):.4g}" if rows else "NA"

    stats = [
        ["features_simulated", len(features)],
        ["features_with_prior", len(priors)],
        [f"features_below_{args.min_reads}_reads", len(features) - len(priors)],
        ["median_H1_prior", median(priors, 1)],
        ["median_H2_prior", median(priors, 2)],
        ["median_H1_misassignment", median(qc, 3)],
        ["median_H2_misassignment", median(qc, 4)],
    ]
    for haplotype, calls in totals.items():
        reads = sum(calls)
        stats.append([f"{haplotype}_simulated_reads", f"{reads:.0f}"])
        for called, value in zip(("H1", "H2", "NonHS"), calls):
            stats.append([f"{haplotype}_reads_called_{called}", f"{value / reads:.4g}" if reads else "NA"])
    write_tsv(args.stats, ["metric", "value"], stats)
    print(f"[build_priors] {len(priors)} features; {len(features) - len(priors)} dropped with fewer than "
          f"{args.min_reads} simulated reads on a haplotype", file=sys.stderr)


if __name__ == "__main__":
    main()
