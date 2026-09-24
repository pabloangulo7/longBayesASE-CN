#!/usr/bin/env python3
"""Classify RNA reads as H1, H2 or NonHS from their alignments and count them.

A read is described by the alignments that reach its best alignment score: when
all of them fall on one haplotype the read is H1 or H2, and when both haplotypes
explain it equally well it is NonHS. The decision rests on the sequence alone,
as in BayesASE, so the simulated reads behind the priors and the real reads are
classified by the same rule whatever the composition of the library.

With --intervals-per-feature, the script also keeps where on its transcript each
counted read lies: a uniform sample of up to that many of the reads counted for
each feature, per level. The mapping priors are simulated from these intervals.

Gene and transcript level are counted in one pass over the BAM. The BAM must keep
every alignment of a read together, as minimap2 writes it; memory then does not
grow with library depth, and a read whose alignments come back in a later block
stops the run instead of being counted twice.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import random
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
LEVELS = ("gene", "transcript")


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
    """Counts, read groups, complex reads and read intervals for one feature level."""

    def __init__(self, prefix: str, intervals_per_feature: int) -> None:
        self.prefix = prefix
        self.counts: dict[str, Counter] = defaultdict(Counter)
        self.complex_reads: dict[str, list[str]] = defaultdict(list)
        self.qc: Counter = Counter()
        self.groups = open(prefix + "_readgroups.tsv", "w")
        self.groups.write("Read_ID\tGroup\n")
        self.intervals_per_feature = intervals_per_feature
        self.intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        self.seen: Counter = Counter()
        # Fixed seed: the same BAM always gives the same intervals.
        self.random = random.Random(1)

    def add(self, read_id: str, h1: list[str], h2: list[str], interval: tuple[str, int, int]) -> None:
        group = classify(h1, h2)
        self.groups.write(f"{read_id}\t{group}\n")
        self.qc[group] += 1
        features = set(h1) | set(h2)
        for feature in features:
            self.counts[feature][group] += 1
            if group in COMPLEX_CATEGORIES:
                self.complex_reads[feature].append(read_id)
        if self.intervals_per_feature and len(features) == 1:
            self.keep_interval(next(iter(features)), interval)

    def keep_interval(self, feature: str, interval: tuple[str, int, int]) -> None:
        """Reservoir sampling: every read of the feature is kept with equal probability."""
        self.seen[feature] += 1
        kept = self.intervals[feature]
        if len(kept) < self.intervals_per_feature:
            kept.append(interval)
        else:
            slot = self.random.randrange(self.seen[feature])
            if slot < self.intervals_per_feature:
                kept[slot] = interval

    def write(self, discarded: Counter) -> None:
        self.groups.close()
        if self.intervals_per_feature:
            with gzip.open(self.prefix + "_read_intervals.tsv.gz", "wt") as out:
                out.write("ID\ttranscript\tstart\tend\n")
                for feature, kept in sorted(self.intervals.items()):
                    for transcript, start, end in sorted(kept):
                        out.write(f"{feature}\t{transcript}\t{start}\t{end}\n")
        with open(self.prefix + "_HS_counts.tsv", "w") as out:
            out.write("ID\t" + "\t".join(CATEGORIES) + "\n")
            # Sorted so that two runs over the same input give the same file.
            for feature, counts in sorted(self.counts.items()):
                out.write(feature + "\t" + "\t".join(str(counts[cat]) for cat in CATEGORIES) + "\n")
        with open(self.prefix + "_complex_reads.tsv", "w") as out:
            out.write("ID\tReads\n")
            for feature, reads in sorted(self.complex_reads.items()):
                out.write(f"{feature}\t{','.join(reads)}\n")
        with open(self.prefix + "_qc.tsv", "w") as out:
            out.write("metric\treads\n")
            for metric, value in sorted((self.qc + discarded).items()):
                out.write(f"{metric}\t{value}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="transcriptome alignments, grouped by read")
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--strand", choices=["fw", "rc", "both"], default="fw",
                        help="orientation a read must have on its transcript")
    parser.add_argument("--intervals-per-feature", type=int, default=0,
                        help="read intervals to keep per feature and level for the priors; 0 keeps none")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    tx2gene = load_tx2gene(args.tx2gene)
    # Sampled per level, so a minor isoform of a highly expressed gene gets its
    # own intervals instead of its share of the gene's.
    counters = {level: LevelCounter(f"{args.output_prefix}.{level}", args.intervals_per_feature)
                for level in LEVELS}
    discarded: Counter = Counter()
    # 8 bytes per read instead of the read names, to check at the end that no
    # read was split into two blocks.
    seen = array("q")

    def finish(read_id: str, best: list[tuple[str, int, int, int]], had_alignment: bool) -> None:
        seen.append(read_key(read_id))
        if not best:
            if had_alignment:
                discarded["Discarded_orientation"] += 1
            return
        h1 = [transcript for transcript, hap, _start, _end in best if hap == 1]
        h2 = [transcript for transcript, hap, _start, _end in best if hap == 2]
        # The read's span on the transcript of its first best alignment.
        interval = (best[0][0], best[0][2], best[0][3])
        counters["transcript"].add(read_id, h1, h2, interval)
        counters["gene"].add(read_id, [gene_of(tx2gene, t) for t in h1], [gene_of(tx2gene, t) for t in h2], interval)

    current = None
    best: list[tuple[str, int, int, int]] = []
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
                best.append((*split_haplotype(alignment.reference_name),
                             alignment.reference_start, alignment.reference_end))
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
