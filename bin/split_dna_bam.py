#!/usr/bin/env python3
"""Split co-optimal DNA alignments, grouped by read, into coverage groups."""

from __future__ import annotations

import argparse
import csv
import itertools
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
    parser.add_argument("--bam", required=True, help="max-score BAM grouped by read")
    parser.add_argument("--readgroups", required=True, help="gene read groups in the same read order")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--hap1-regex", default=r"_hap1(?:$|_)")
    parser.add_argument("--hap2-regex", default=r"_hap2(?:$|_)")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    hap1_pattern = re.compile(args.hap1_regex)
    hap2_pattern = re.compile(args.hap2_regex)

    unclassified = 0
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as source, open(args.readgroups) as groups:
        readgroups = csv.DictReader(groups, delimiter="\t")
        outputs = {
            name: pysam.AlignmentFile(f"{args.output_prefix}.{name}.bam", "wb", template=source, threads=args.threads)
            for name in ("HS", "NonHS", "HS_copy", "NonHS_copy")
        }
        try:
            for read_name, alignments in itertools.groupby(source, key=lambda alignment: alignment.query_name):
                row = next(readgroups, None)
                if row is None or row["Read_ID"] != read_name:
                    raise ValueError(f"read {read_name}: the BAM and the read-group table are not in the same read order")
                destination = target_group(row["Group"])
                if destination is None:
                    continue
                # One alignment per haplotype: the first in coordinate order.
                kept: dict[str, pysam.AlignedSegment] = {}
                for alignment in alignments:
                    reference = source.get_reference_name(alignment.reference_id)
                    if hap1_pattern.search(reference):
                        hap = "H1"
                    elif hap2_pattern.search(reference):
                        hap = "H2"
                    else:
                        unclassified += 1
                        continue
                    position = (alignment.reference_id, alignment.reference_start)
                    if hap not in kept or position < (kept[hap].reference_id, kept[hap].reference_start):
                        kept[hap] = alignment
                for alignment in kept.values():
                    outputs[destination].write(alignment)
        finally:
            for handle in outputs.values():
                handle.close()
    print(f"[split_dna_bam] dropped {unclassified} alignments on unclassified contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
