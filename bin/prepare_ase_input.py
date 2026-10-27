#!/usr/bin/env python3
"""Join haplotype-specific counts, copy number and mapping priors into the ASE input.

One row per library and feature, with the library's condition as its group.
The misassignment rates of the priors are carried along so that the diffASE
step can leave out the features they make unreliable.
"""

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


def load_priors(path: str) -> dict[str, tuple[float, float, float, float]]:
    """r1, r2 and the two misassignment rates of every feature (0 when the table has none)."""
    priors = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not {"ID", "H1_prior", "H2_prior"}.issubset(fields):
            raise ValueError("priors require columns: ID, H1_prior, H2_prior")
        has_misassignment = {"H1_misassignment", "H2_misassignment"}.issubset(fields)
        for row in reader:
            r1, r2 = float(row["H1_prior"]), float(row["H2_prior"])
            if not 0 < r1 < 1 or not 0 < r2 < 1:
                raise ValueError(f"feature {row['ID']}: priors must be strictly between 0 and 1")
            wrong = ((float(row["H1_misassignment"]), float(row["H2_misassignment"]))
                     if has_misassignment else (0.0, 0.0))
            priors[row["ID"]] = r1, r2, *wrong
    return priors


def load_copy_number(path: str) -> dict[tuple[str, str], tuple[float, float]]:
    result = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample", "ID", "CN_H1", "CN_H2"}.issubset(reader.fieldnames or []):
            raise ValueError("copy number requires columns: sample, ID, CN_H1, CN_H2")
        for row in reader:
            values = float(row["CN_H1"]), float(row["CN_H2"])
            if not all(math.isfinite(value) and value >= 0 for value in values):
                raise ValueError(f"gene {row['ID']}: copy number must be finite and >= 0")
            result[(row["sample"], row["ID"])] = values
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--hs-counts", required=True)
    parser.add_argument("--copy-number", required=True)
    parser.add_argument("--priors", required=True)
    parser.add_argument("--tx2gene", help="maps transcript-level counts to the gene their copy number is kept under")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    samples = load_samples(args.samplesheet)
    priors = load_priors(args.priors)
    copy_number = load_copy_number(args.copy_number)
    tx2gene = load_tx2gene(args.tx2gene)

    def gene_copy_number(sample: str, feature: str) -> tuple[float, float] | None:
        # Copy number is keyed by the DNA sample and the gene.
        return copy_number.get((samples[sample]["dna_id"], tx2gene.get(feature, feature)))

    def model_row(sample: str, feature: str, counts: tuple[int, int, int],
                  cn: tuple[float, float]) -> dict[str, object]:
        r1, r2, wrong1, wrong2 = priors[feature]
        return {
            "sample": sample, "dna_id": samples[sample]["dna_id"], "group": samples[sample]["condition"],
            "ID": feature, "H1_counts": counts[0], "H2_counts": counts[1], "NonHS_counts": counts[2],
            "CN_H1": cn[0], "CN_H2": cn[1], "H1_prior": r1, "H2_prior": r2,
            "H1_misassignment": wrong1, "H2_misassignment": wrong2,
        }

    rows_out = []
    missing_cn, missing_prior = set(), set()
    libraries: set[str] = set()
    counted: set[tuple[str, str]] = set()
    with open(args.hs_counts) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample", "ID", "H1", "H2", "NonHS"}.issubset(reader.fieldnames or []):
            raise ValueError("HS counts require columns: sample, ID, H1, H2, NonHS")
        for row in reader:
            sample, feature = row["sample"], row["ID"]
            if sample not in samples:
                raise ValueError(f"HS counts have no samplesheet row: {sample}")
            libraries.add(sample)
            cn = gene_copy_number(sample, feature)
            if cn is None:
                missing_cn.add(feature)
                continue
            if feature not in priors:
                missing_prior.add(feature)
                continue
            counted.add((sample, feature))
            counts = number(row, "H1"), number(row, "H2"), number(row, "NonHS")
            rows_out.append(model_row(sample, feature, counts, cn))

    # A library with no read for a feature has no row in the count table, but it
    # did measure zero; without the row the model would drop that library.
    filled = 0
    for feature in sorted({feature for _sample, feature in counted}):
        for sample in sorted(libraries):
            if (sample, feature) in counted:
                continue
            cn = gene_copy_number(sample, feature)
            if cn is not None:
                rows_out.append(model_row(sample, feature, (0, 0, 0), cn))
                filled += 1
    print(
        f"[prepare_ase_input] dropped {len(missing_cn | missing_prior)} features: "
        f"{len(missing_cn)} without copy number, {len(missing_prior)} without prior; "
        f"added {filled} zero-count rows",
        file=sys.stderr,
    )
    if not rows_out:
        raise ValueError("no features remained after joining counts, copy number and priors")
    rows_out.sort(key=lambda row: (row["ID"], row["group"], row["sample"]))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
