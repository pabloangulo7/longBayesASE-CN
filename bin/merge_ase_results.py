#!/usr/bin/env python3
"""Merge the ASE shard results and flag the features the priors make unreliable.

Every feature is fitted; filter_flag marks those whose result should not be
trusted: low_prior when min(r1, r2) <= --min-prior (few of its reads cover a
heterozygous site on one haplotype, so r1/r2 and the bias correction are
poorly determined), and misassignment when simulated reads come out on the
wrong haplotype more often than --max-misassignment, which the model has no
term for. pass otherwise.
"""

from __future__ import annotations

import argparse
import csv


def load_flags(path: str, min_prior: float, max_misassignment: float) -> dict[str, str]:
    flags = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["ID"] in flags:
                continue
            reasons = []
            if min(float(row["H1_prior"]), float(row["H2_prior"])) <= min_prior:
                reasons.append("low_prior")
            wrong = max(float(row.get("H1_misassignment") or 0), float(row.get("H2_misassignment") or 0))
            if wrong > max_misassignment:
                reasons.append("misassignment")
            flags[row["ID"]] = ";".join(reasons) or "pass"
    return flags


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--ase-input", required=True)
    parser.add_argument("--min-prior", type=float, default=0.01)
    parser.add_argument("--max-misassignment", type=float, default=0.05)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    flags = load_flags(args.ase_input, args.min_prior, args.max_misassignment)
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
    for row in rows:
        row["filter_flag"] = flags[row["ID"]]
    rows.sort(key=lambda row: (row.get("test", ""), row.get("ID", "")))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*(fields or ["ID", "test"]), "filter_flag"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
