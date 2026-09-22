#!/usr/bin/env python3
"""Join unified RNA counts, copy number, priors and sample metadata."""

from __future__ import annotations

import argparse
import csv
import math
import sys

from common import load_tx2gene


def number(row: dict[str, str], name: str) -> int:
    return int(float(row.get(name) or 0))


def load_samples(path: str) -> dict[str, dict[str, str]]:
    with open(path) as handle:
        return {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}


def uncommented(handle):
    """Drop the provenance header written by build_priors.py."""
    return (line for line in handle if not line.startswith("#"))


def load_priors(path: str) -> dict[str, tuple[float, float]]:
    priors = {}
    with open(path) as handle:
        reader = csv.DictReader(uncommented(handle), delimiter="\t")
        if not {"ID", "H1_prior", "H2_prior"}.issubset(reader.fieldnames or []):
            raise ValueError("priors require columns: ID, H1_prior, H2_prior")
        for row in reader:
            values = float(row["H1_prior"]), float(row["H2_prior"])
            if not all(0 < value < 1 for value in values):
                raise ValueError(f"gene {row['ID']}: priors must be strictly between 0 and 1")
            priors[row["ID"]] = values
    return priors


def load_copy_number(path: str) -> dict[tuple[str, str], tuple[float, float]]:
    result = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        sample_column = "sample" if "sample" in fields else "dna_id" if "dna_id" in fields else ""
        if not sample_column or not {"ID", "CN_H1", "CN_H2"}.issubset(fields):
            raise ValueError("copy number requires columns: sample, ID, CN_H1, CN_H2")
        for row in reader:
            values = float(row["CN_H1"]), float(row["CN_H2"])
            if not all(math.isfinite(value) and value >= 0 for value in values):
                raise ValueError(f"gene {row['ID']}: copy number must be finite and >= 0")
            result[(row[sample_column], row["ID"])] = values
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--rna-counts")
    group.add_argument("--counts")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--copy-number")
    parser.add_argument("--priors", required=True)
    parser.add_argument("--tx2gene", help="used when the counts are isoform-level")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.rna_counts and not args.copy_number:
        parser.error("--rna-counts requires --copy-number; use --counts for a combined table")

    samples = load_samples(args.samplesheet)
    priors = load_priors(args.priors)
    copy_number = load_copy_number(args.copy_number) if args.copy_number else {}
    tx2gene = load_tx2gene(args.tx2gene)
    rows_out = []
    missing_cn, missing_prior = set(), set()
    with open(args.counts or args.rna_counts) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not {"sample", "ID", "H1", "H2", "NonHS"}.issubset(fields):
            raise ValueError("RNA counts require columns: sample, ID, H1, H2, NonHS")
        if args.counts and not {"CN_H1", "CN_H2"}.issubset(fields):
            raise ValueError("combined counts also require columns: CN_H1, CN_H2")
        for row in reader:
            sample, gene = row["sample"], row["ID"]
            if sample not in samples:
                raise ValueError(f"RNA counts have no samplesheet row: {sample}")
            meta = samples[sample]
            # Copy number is always keyed by gene; isoform-level counts reach it
            # through tx2gene.
            cn = ((float(row["CN_H1"]), float(row["CN_H2"])) if args.counts
                  else copy_number.get((meta["dna_id"], tx2gene.get(gene, gene))))
            prior = priors.get(gene)
            if cn is None:
                missing_cn.add(gene)
                continue
            if prior is None:
                missing_prior.add(gene)
                continue
            rows_out.append({
                "sample": sample, "dna_id": meta["dna_id"], "group": meta["condition"],
                "ID": gene,
                "H1_counts": number(row, "H1"), "H2_counts": number(row, "H2"),
                "NonHS_counts": number(row, "NonHS"),
                "CN_H1": cn[0], "CN_H2": cn[1],
                "H1_prior": prior[0], "H2_prior": prior[1],
            })
    print(
        f"[prepare_ase_input] dropped {len(missing_cn | missing_prior)} genes: "
        f"{len(missing_cn)} missing copy number, {len(missing_prior)} missing prior",
        file=sys.stderr,
    )
    if not rows_out:
        raise ValueError("no genes remained after joining counts, copy number and priors")
    rows_out.sort(key=lambda row: (row["ID"], row["group"], row["sample"]))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
