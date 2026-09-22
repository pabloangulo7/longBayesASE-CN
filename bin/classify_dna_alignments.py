#!/usr/bin/env python3
"""Keep co-optimal DNA alignments and classify their haplotype behaviour."""

from __future__ import annotations

import argparse
import itertools
import re

import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="query-name collated BAM")
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--readgroups", required=True)
    parser.add_argument("--hap1-regex", default=r"_hap1(?:$|_)")
    parser.add_argument("--hap2-regex", default=r"_hap2(?:$|_)")
    parser.add_argument("--score-tag", default="AS")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    hap1_pattern = re.compile(args.hap1_regex)
    hap2_pattern = re.compile(args.hap2_regex)

    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as source, pysam.AlignmentFile(
        args.output_bam, "wb", template=source, threads=args.threads
    ) as target, open(args.readgroups, "w") as table:
        table.write("Read_ID\tH1\tH2\tGroup\tMax_Score\tUnclassified\n")
        for read_name, alignments_iter in itertools.groupby(source, key=lambda item: item.query_name):
            alignments = [
                aln
                for aln in alignments_iter
                if not aln.is_unmapped and not aln.is_supplementary and aln.has_tag(args.score_tag)
            ]
            if not alignments:
                continue
            max_score = max(aln.get_tag(args.score_tag) for aln in alignments)
            best = [aln for aln in alignments if aln.get_tag(args.score_tag) == max_score]
            hap_counts = {"H1": 0, "H2": 0, "Unclassified": 0}
            for aln in best:
                reference = source.get_reference_name(aln.reference_id)
                if hap1_pattern.search(reference):
                    hap = "H1"
                elif hap2_pattern.search(reference):
                    hap = "H2"
                else:
                    hap = "Unclassified"
                hap_counts[hap] += 1
                target.write(aln)
            if hap_counts["H1"] and hap_counts["H2"]:
                group = "NonHS"
            elif hap_counts["H1"]:
                group = "H1"
            elif hap_counts["H2"]:
                group = "H2"
            else:
                group = "Unclassified"
            if hap_counts["H1"] > 1 or hap_counts["H2"] > 1:
                group += "_multimapping"
            table.write(
                f"{read_name}\t{hap_counts['H1']}\t{hap_counts['H2']}\t{group}\t{max_score}"
                f"\t{hap_counts['Unclassified']}\n"
            )


if __name__ == "__main__":
    main()
