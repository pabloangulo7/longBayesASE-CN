#!/usr/bin/env python3
"""Create exhaustive deterministic FASTQ tiles from every transcript.

Without --isoform-usage every transcript is tiled once, so the per-gene
assignment probability is an average over isoforms weighted by their length.
With --isoform-usage each tile is emitted k times, with k proportional to how
much that isoform is used relative to the dominant isoform of its own gene, so
the average is weighted by expression instead. Only within-gene ratios matter,
because the priors are estimated per gene.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict

from common import HAP_SUFFIX_RE, fasta_records, load_tx2gene, open_text


def tile_starts(length: int, read_length: int, step: int) -> list[int]:
    if length <= read_length:
        return [0]
    last = length - read_length
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    return starts


def load_usage(path: str | None) -> dict[str, float]:
    """Sum the Oarfish transcript table (ID, num_reads) over samples.

    Summing here is deliberate: the prior has to be identical in every
    condition, so it can never be built from one condition's expression.
    """
    weights: dict[str, float] = defaultdict(float)
    if not path:
        return weights
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not fields:
            return weights
        if not {"ID", "num_reads"}.issubset(fields):
            raise ValueError("isoform usage table requires columns: ID, num_reads")
        for row in reader:
            weights[row["ID"]] += float(row["num_reads"])
    return weights


def replicate_counts(
    transcripts: list[str],
    weights: dict[str, float],
    tx2gene: dict[str, str],
    max_replicates: int,
) -> dict[str, int]:
    """Replicates per transcript, normalised against its own gene's maximum."""
    if not weights or max_replicates <= 1:
        return {transcript: 1 for transcript in transcripts}
    bases = {transcript: HAP_SUFFIX_RE.sub("", transcript) for transcript in transcripts}
    top: dict[str, float] = defaultdict(float)
    for base in bases.values():
        gene = tx2gene.get(base, base)
        top[gene] = max(top[gene], weights.get(base, 0.0))
    counts = {}
    for transcript, base in bases.items():
        gene_top = top[tx2gene.get(base, base)]
        share = weights.get(base, 0.0) / gene_top if gene_top > 0 else 0.0
        counts[transcript] = min(max(1, round(max_replicates * share)), max_replicates)
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcriptome", required=True)
    parser.add_argument("--read-length", type=int, default=1000)
    parser.add_argument("--step", type=int, default=100)
    parser.add_argument("--isoform-usage", help="Oarfish transcript table with columns ID and num_reads")
    parser.add_argument("--tx2gene", help="transcript_id/gene_id table used to group isoforms")
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--output", required=True)
    parser.add_argument("--weighting-report", help="write the applied weighting mode here")
    args = parser.parse_args()
    if args.read_length < 1 or args.step < 1:
        parser.error("--read-length and --step must be positive")
    if args.max_replicates < 1:
        parser.error("--max-replicates must be at least 1")

    weights = load_usage(args.isoform_usage)
    tx2gene = load_tx2gene(args.tx2gene)
    # A first pass reads only the headers, so the sequences are never all in memory.
    with open_text(args.transcriptome) as handle:
        identifiers = [line[1:].split()[0] for line in handle if line.startswith(">")]
    replicates = replicate_counts(identifiers, weights, tx2gene, args.max_replicates)
    weighted = bool(weights) and args.max_replicates > 1
    mode = (f"expression (max_replicates={args.max_replicates}, "
            f"{len(weights)} transcripts with usage)" if weighted else "uniform")

    opener = gzip.open if args.output.endswith(".gz") else open
    written = 0
    with opener(args.output, "wt") as output:
        for transcript, _header, sequence in fasta_records(args.transcriptome):
            sequence = sequence.upper()
            if not sequence:
                continue
            safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", transcript)
            for start in tile_starts(len(sequence), args.read_length, args.step):
                tile = sequence[start : start + args.read_length]
                name = f"{safe_id}__tile_{start + 1}_{start + len(tile)}"
                for replicate in range(replicates.get(transcript, 1)):
                    suffix = "" if replicate == 0 else f"__rep{replicate + 1}"
                    output.write(f"@{name}{suffix}\n{tile}\n+\n{'I' * len(tile)}\n")
                    written += 1
    if written == 0:
        raise ValueError("transcriptome contains no non-empty sequences")
    if args.weighting_report:
        with open(args.weighting_report, "w") as handle:
            handle.write(mode + "\n")
    print(f"[tile_transcripts] {written} tiles from {len(identifiers)} transcripts | weighting: {mode}")


if __name__ == "__main__":
    main()
