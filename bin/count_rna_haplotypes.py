#!/usr/bin/env python3
"""Classify RNA reads as H1, H2 or NonHS from their alignments and count them.

A read is described by the alignments that reach its best alignment score: when
all of them fall on one haplotype the read is H1 or H2, and when both haplotypes
explain it equally well it is NonHS. The decision rests on the sequence alone,
as in BayesASE, so the simulated reads behind the priors and the real reads are
classified by the same rule whatever the composition of the library.

Gene and isoform level are counted in one pass over the BAM. The BAM must keep
every alignment of a read together, as minimap2 writes it; memory then does not
grow with library depth, and a read whose alignments come back in a later block
stops the run instead of being counted twice.
"""

from __future__ import annotations

import argparse
import hashlib
from array import array
from collections import Counter, defaultdict

import numpy as np
import pysam

from common import gene_of, load_tx2gene, split_haplotype


CATEGORIES = [
    "H1", "H1_multimapping", "H1_multimapping_multigene",
    "H2", "H2_multimapping", "H2_multimapping_multigene",
    "NonHS", "NonHS_multimapping", "NonHS_multimapping_multigene",
    "NonHS_complex", "NonHS_multimapping_complex",
]
COMPLEX_CATEGORIES = {
    "H1_multimapping_multigene", "H2_multimapping_multigene",
    "NonHS_multimapping_multigene", "NonHS_complex", "NonHS_multimapping_complex",
}
LEVELS = ("genes", "isoforms")


def read_key(read_id: str) -> int:
    return int.from_bytes(hashlib.blake2b(read_id.encode(), digest_size=8).digest(), "little", signed=True)


def classify(h1: list[str], h2: list[str]) -> str:
    """Name a read's group from the features of its best alignments on each haplotype."""
    if h1 and h2:
        multimapping = len(h1) > 1 or len(h2) > 1
        # Both haplotypes explain the read, but not through the same features.
        if set(h1) != set(h2):
            return "NonHS_multimapping_complex" if multimapping else "NonHS_complex"
        haplotype = "NonHS"
    else:
        haplotype = "H1" if h1 else "H2"
        multimapping = len(h1) + len(h2) > 1
    if not multimapping:
        return haplotype
    if len(set(h1) | set(h2)) > 1:
        return f"{haplotype}_multimapping_multigene"
    return f"{haplotype}_multimapping"


class LevelCounter:
    """Counts, read groups and complex reads for one feature level."""

    def __init__(self, prefix: str) -> None:
        self.prefix = prefix
        self.counts: dict[str, Counter] = defaultdict(Counter)
        self.complex_reads: dict[str, list[str]] = defaultdict(list)
        self.qc: Counter = Counter()
        self.groups = open(prefix + ".readgroups.tsv", "w")
        self.groups.write("Read_ID\tGroup\n")

    def add(self, read_id: str, h1: list[str], h2: list[str]) -> None:
        group = classify(h1, h2)
        self.groups.write(f"{read_id}\t{group}\n")
        self.qc[group] += 1
        for feature in set(h1) | set(h2):
            self.counts[feature][group] += 1
            if group in COMPLEX_CATEGORIES:
                self.complex_reads[feature].append(read_id)

    def write(self, discarded: Counter) -> None:
        self.groups.close()
        with open(self.prefix + ".counts.tsv", "w") as out:
            out.write("ID\t" + "\t".join(CATEGORIES) + "\n")
            # Sorted so that two runs over the same input give the same file.
            for feature, counts in sorted(self.counts.items()):
                out.write(feature + "\t" + "\t".join(str(counts[cat]) for cat in CATEGORIES) + "\n")
        with open(self.prefix + ".complex_reads.tsv", "w") as out:
            out.write("ID\tReads\n")
            for feature, reads in sorted(self.complex_reads.items()):
                out.write(f"{feature}\t{','.join(reads)}\n")
        with open(self.prefix + ".qc.tsv", "w") as out:
            out.write("metric\treads\n")
            for metric, value in sorted((self.qc + discarded).items()):
                out.write(f"{metric}\t{value}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="transcriptome alignments, grouped by read")
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--strand", choices=["fw", "rc", "both"], default="fw",
                        help="orientation a read must have on its transcript")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    tx2gene = load_tx2gene(args.tx2gene)
    counters = {level: LevelCounter(f"{args.output_prefix}.{level}") for level in LEVELS}
    discarded: Counter = Counter()
    # 8 bytes per read instead of the read names, to check at the end that no
    # read was split into two blocks.
    seen = array("q")

    def finish(read_id: str, best: list[tuple[str, int]], had_alignment: bool) -> None:
        seen.append(read_key(read_id))
        if not best:
            if had_alignment:
                discarded["Discarded_orientation"] += 1
            return
        h1 = [transcript for transcript, hap in best if hap == 1]
        h2 = [transcript for transcript, hap in best if hap == 2]
        counters["isoforms"].add(read_id, h1, h2)
        counters["genes"].add(read_id, [gene_of(tx2gene, t) for t in h1], [gene_of(tx2gene, t) for t in h2])

    current = None
    best: list[tuple[str, int]] = []
    best_score = None
    had_alignment = False
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False, threads=args.threads) as bam:
        for alignment in bam.fetch(until_eof=True):
            read_id = alignment.query_name
            if read_id != current:
                if current is not None:
                    finish(current, best, had_alignment)
                current, best, best_score, had_alignment = read_id, [], None, False
            if alignment.is_unmapped or alignment.is_supplementary:
                continue
            had_alignment = True
            if args.strand != "both" and alignment.is_reverse != (args.strand == "rc"):
                continue
            score = alignment.get_tag("AS")
            if best_score is None or score > best_score:
                best, best_score = [], score
            if score == best_score:
                best.append(split_haplotype(alignment.reference_name))
        if current is not None:
            finish(current, best, had_alignment)

    keys = np.frombuffer(seen, dtype=np.int64)
    keys.sort()
    if keys.size > 1 and np.any(keys[1:] == keys[:-1]):
        raise ValueError(f"{args.bam}: the alignments of at least one read are not together; "
                         "use the aligner's output order, not a coordinate-sorted BAM")
    for counter in counters.values():
        counter.write(discarded)


if __name__ == "__main__":
    main()
