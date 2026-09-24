#!/usr/bin/env python3
"""Simulate the reads behind the mapping priors from where the real reads lie.

Every read interval sampled from the RNA libraries (transcript, start, end) is cut
from one haplotype's copy of that transcript, so both haplotypes are simulated at
exactly the positions the libraries cover. Whether a read overlaps a
heterozygous site depends on where it starts and ends, so these positions, and
not the reads' alleles, are what the priors need from the data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from array import array
from collections import defaultdict

from common import fasta_records, open_text, split_haplotype

# Intervals shorter than this are not worth aligning.
MIN_LENGTH = 20


def load_intervals(path: str) -> dict[str, array]:
    """Start and end of every interval, flattened per transcript."""
    by_transcript: dict[str, array] = defaultdict(lambda: array("i"))
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"transcript", "start", "end"}.issubset(reader.fieldnames or []):
            raise ValueError("read intervals require columns: transcript, start, end")
        for row in reader:
            by_transcript[row["transcript"]].extend((int(row["start"]), int(row["end"])))
    return by_transcript


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--intervals", required=True, help="read intervals from the RNA libraries")
    parser.add_argument("--transcriptome", required=True, help="transcripts of one haplotype, named <transcript>_hap1|2")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    intervals = load_intervals(args.intervals)
    if not intervals:
        raise ValueError(f"{args.intervals} holds no read intervals")
    opener = gzip.open if args.output.endswith(".gz") else open
    written = skipped = 0
    with opener(args.output, "wt") as output:
        for name, _header, sequence in fasta_records(args.transcriptome):
            spans = intervals.get(split_haplotype(name)[0])
            if not spans:
                continue
            sequence = sequence.upper()
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            # Coordinates come from whichever haplotype the real read aligned to,
            # so they are clipped to this copy's length.
            for index in range(0, len(spans), 2):
                start, end = spans[index], min(spans[index + 1], len(sequence))
                if end - start < MIN_LENGTH:
                    skipped += 1
                    continue
                read = sequence[start:end]
                output.write(f"@{safe}__{start + 1}_{end}__{index // 2}\n{read}\n+\n{'I' * len(read)}\n")
                written += 1
    if written == 0:
        raise ValueError(f"no read interval matched a transcript of {args.transcriptome}")
    print(f"[simulate_prior_reads] {written} reads from {len(intervals)} transcripts; "
          f"{skipped} intervals under {MIN_LENGTH} bp skipped")


if __name__ == "__main__":
    main()
