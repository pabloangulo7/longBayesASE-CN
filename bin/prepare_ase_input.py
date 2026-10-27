#!/usr/bin/env python3
"""Join haplotype-specific counts, copy number, priors and sample metadata.

Genes whose mapping priors cannot inform the model are left out here: those
whose reads rarely cover a heterozygous site on one of the haplotypes
(min(r1, r2) <= --min-prior), and those whose simulated reads come out on the
wrong haplotype more often than --max-misassignment, which the model has no
term for.
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


def load_priors(path: str) -> dict[str, tuple[float, float, float]]:
    """r1, r2 and the larger of the two misassignment rates of every gene.

    Priors supplied without misassignment columns are taken to have none.
    """
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
                raise ValueError(f"gene {row['ID']}: priors must be strictly between 0 and 1")
            misassignment = (max(float(row["H1_misassignment"]), float(row["H2_misassignment"]))
                             if has_misassignment else 0.0)
            priors[row["ID"]] = r1, r2, misassignment
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
    group.add_argument("--hs-counts")
    group.add_argument("--counts")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--copy-number")
    parser.add_argument("--priors", required=True)
    parser.add_argument("--tx2gene", help="used when the counts are transcript-level")
    parser.add_argument("--min-prior", type=float, default=0.01,
                        help="genes with min(r1, r2) at or below this are left out")
    parser.add_argument("--max-misassignment", type=float, default=0.05,
                        help="genes whose simulated reads land on the wrong haplotype more often are left out")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.hs_counts and not args.copy_number:
        parser.error("--hs-counts requires --copy-number; use --counts for a combined table")

    samples = load_samples(args.samplesheet)
    priors = load_priors(args.priors)
    uninformative = {gene for gene, (r1, r2, _) in priors.items() if min(r1, r2) <= args.min_prior}
    misassigned = {gene for gene, (_, _, rate) in priors.items()
                   if rate > args.max_misassignment} - uninformative
    copy_number = load_copy_number(args.copy_number) if args.copy_number else {}
    tx2gene = load_tx2gene(args.tx2gene)

    def model_row(sample: str, gene: str, counts: tuple[int, int, int],
                  cn: tuple[float, float]) -> dict[str, object]:
        meta = samples[sample]
        return {
            "sample": sample, "dna_id": meta["dna_id"], "group": meta["condition"], "ID": gene,
            "H1_counts": counts[0], "H2_counts": counts[1], "NonHS_counts": counts[2],
            "CN_H1": cn[0], "CN_H2": cn[1], "H1_prior": priors[gene][0], "H2_prior": priors[gene][1],
        }

    def table_copy_number(sample: str, gene: str) -> tuple[float, float] | None:
        # Copy number is always keyed by gene; transcript-level counts reach it
        # through tx2gene.
        return copy_number.get((samples[sample]["dna_id"], tx2gene.get(gene, gene)))

    rows_out = []
    missing_cn, missing_prior, filtered = set(), set(), set()
    libraries: set[str] = set()
    kept_cn: dict[tuple[str, str], tuple[float, float]] = {}
    with open(args.counts or args.hs_counts) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not {"sample", "ID", "H1", "H2", "NonHS"}.issubset(fields):
            raise ValueError("HS counts require columns: sample, ID, H1, H2, NonHS")
        if args.counts and not {"CN_H1", "CN_H2"}.issubset(fields):
            raise ValueError("combined counts also require columns: CN_H1, CN_H2")
        for row in reader:
            sample, gene = row["sample"], row["ID"]
            if sample not in samples:
                raise ValueError(f"HS counts have no samplesheet row: {sample}")
            libraries.add(sample)
            if gene in uninformative or gene in misassigned:
                filtered.add(gene)
                continue
            cn = ((float(row["CN_H1"]), float(row["CN_H2"])) if args.counts
                  else table_copy_number(sample, gene))
            if cn is None:
                missing_cn.add(gene)
                continue
            if gene not in priors:
                missing_prior.add(gene)
                continue
            kept_cn[(sample, gene)] = cn
            counts = number(row, "H1"), number(row, "H2"), number(row, "NonHS")
            rows_out.append(model_row(sample, gene, counts, cn))

    # A library with no read for a gene has no row in the count table, but it
    # did measure zero; without the row the model would drop that library. A
    # combined table carries copy number only on its rows, so there it comes
    # from another library of the same DNA sample.
    cn_by_dna = {(samples[sample]["dna_id"], gene): cn for (sample, gene), cn in kept_cn.items()}
    filled = unfillable = 0
    for gene in sorted({gene for _sample, gene in kept_cn}):
        for sample in sorted(libraries):
            if (sample, gene) in kept_cn:
                continue
            cn = (cn_by_dna.get((samples[sample]["dna_id"], gene)) if args.counts
                  else table_copy_number(sample, gene))
            if cn is None:
                unfillable += 1
                continue
            rows_out.append(model_row(sample, gene, (0, 0, 0), cn))
            filled += 1
    print(
        f"[prepare_ase_input] dropped {len(missing_cn | missing_prior)} genes: "
        f"{len(missing_cn)} missing copy number, {len(missing_prior)} missing prior; "
        f"left out {len(filtered & uninformative)} with min(r1, r2) <= {args.min_prior} and "
        f"{len(filtered & misassigned)} with misassignment > {args.max_misassignment}; "
        f"added {filled} zero-count rows, {unfillable} left out for lack of copy number",
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
