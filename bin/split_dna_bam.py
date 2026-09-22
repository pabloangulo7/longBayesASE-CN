#!/usr/bin/env python3
"""Split coordinate-sorted co-optimal DNA alignments into coverage groups."""

from __future__ import annotations

import argparse
import csv
import re
import sys

import pysam


VALID_READGROUPS = {
    "H1", "H2", "NonHS",
    "H1_multigene", "H2_multigene", "NonHS_multigene",
    "H1_multimapping", "H2_multimapping", "NonHS_multimapping",
    "H1_multimapping_multigene", "H2_multimapping_multigene", "NonHS_multimapping_multigene",
    "H1_multimapping_multigene_copy", "H2_multimapping_multigene_copy",
    "NonHS_multigene_copy", "NonHS_multimapping_multigene_copy",
}


def read_groups(path: str) -> dict[str, str]:
    with open(path) as handle:
        return {
            row["Read_ID"]: row["Group"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def target_group(group: str) -> str | None:
    if group not in VALID_READGROUPS:
        return None
    is_hs = group.startswith("H1") or group.startswith("H2")
    is_copy = "_copy" in group
    if is_hs:
        return "HS_copy" if is_copy else "HS"
    if group.startswith("NonHS"):
        return "NonHS_copy" if is_copy else "NonHS"
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="coordinate-sorted max-score BAM")
    parser.add_argument("--readgroups", required=True, help="gene read-group TSV")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--hap1-regex", default=r"_hap1(?:$|_)")
    parser.add_argument("--hap2-regex", default=r"_hap2(?:$|_)")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    hap1_pattern = re.compile(args.hap1_regex)
    hap2_pattern = re.compile(args.hap2_regex)

    destinations = {
        read_id: destination
        for read_id, group in read_groups(args.readgroups).items()
        if (destination := target_group(group)) is not None
    }
    unclassified = 0
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as source:
        outputs = {
            name: pysam.AlignmentFile(f"{args.output_prefix}.{name}.bam", "wb", template=source, threads=args.threads)
            for name in ("HS", "NonHS", "HS_copy", "NonHS_copy")
        }
        try:
            seen_haplotypes: dict[str, set[str]] = {"H1": set(), "H2": set()}
            for alignment in source:
                read_name = alignment.query_name
                destination = destinations.get(read_name)
                if destination is None:
                    continue
                reference = source.get_reference_name(alignment.reference_id)
                if hap1_pattern.search(reference):
                    hap = "H1"
                elif hap2_pattern.search(reference):
                    hap = "H2"
                else:
                    unclassified += 1
                    continue
                if read_name in seen_haplotypes[hap]:
                    continue
                outputs[destination].write(alignment)
                seen_haplotypes[hap].add(read_name)
        finally:
            for handle in outputs.values():
                handle.close()
    print(f"[split_dna_bam] dropped {unclassified} alignments on unclassified contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
