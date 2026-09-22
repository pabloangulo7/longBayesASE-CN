#!/usr/bin/env python3
"""Estimate per-gene haplotype assignment probabilities from simulations."""

from __future__ import annotations

import argparse
import csv
import sys


def read_counts(path: str) -> dict[str, dict[str, float]]:
    output = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values = {name: float(row.get(name) or 0) for name in row if name != "ID"}
            output[row["ID"]] = values
    return output


def collapsed(row: dict[str, float]) -> tuple[float, float, float]:
    return (
        row.get("H1", 0) + row.get("H1_multimapping", 0),
        row.get("H2", 0) + row.get("H2_multimapping", 0),
        row.get("NonHS", 0) + row.get("NonHS_multimapping", 0),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hap1-counts", required=True)
    parser.add_argument("--hap2-counts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-simulated-reads", type=int, default=1)
    parser.add_argument("--epsilon", type=float, default=0.001)
    parser.add_argument("--weighting", default="uniform",
                        help="how the tiles were weighted; recorded in the output header")
    args = parser.parse_args()
    h1_rows, h2_rows = read_counts(args.hap1_counts), read_counts(args.hap2_counts)
    shared = sorted(set(h1_rows) & set(h2_rows))
    dropped = 0
    with open(args.output, "w") as handle:
        # Provenance: the same data give different priors under different
        # weighting, so the file has to say which one produced it.
        handle.write(f"# weighting: {args.weighting}\n")
        handle.write("ID\tH1_prior\tH2_prior\tH1_misassignment\tH2_misassignment\tn_hap1\tn_hap2\n")
        for gene in shared:
            h1, h1_wrong, h1_nonhs = collapsed(h1_rows[gene])
            h2_wrong, h2, h2_nonhs = collapsed(h2_rows[gene])
            n1, n2 = h1 + h1_wrong + h1_nonhs, h2 + h2_wrong + h2_nonhs
            if n1 < args.min_simulated_reads or n2 < args.min_simulated_reads:
                dropped += 1
                continue
            prior1 = min(max(h1 / n1, args.epsilon), 1 - args.epsilon)
            prior2 = min(max(h2 / n2, args.epsilon), 1 - args.epsilon)
            handle.write(
                f"{gene}\t{prior1:.8g}\t{prior2:.8g}\t{h1_wrong/n1:.8g}\t{h2_wrong/n2:.8g}\t{n1:.0f}\t{n2:.0f}\n"
            )
    print(f"[build_priors] dropped {dropped} genes below --min-simulated-reads", file=sys.stderr)


if __name__ == "__main__":
    main()
