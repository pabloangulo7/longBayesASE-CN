#!/usr/bin/env python3
"""Split a _hap1/_hap2 diploid assembly and build Liftoff chromosome maps."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from common import fasta_records, write_fasta_record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--hap1-regex", default=r"_hap1$")
    parser.add_argument("--hap2-regex", default=r"_hap2$")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    patterns = {"hap1": re.compile(args.hap1_regex), "hap2": re.compile(args.hap2_regex)}
    counts = {"hap1": 0, "hap2": 0}
    fasta_handles = {hap: open(output / f"{hap}.fa", "w") for hap in patterns}
    map_handles = {hap: open(output / f"{hap}.chroms.csv", "w") for hap in patterns}
    excluded = open(output / "excluded_contigs.tsv", "w")
    excluded.write("contig\treason\n")
    try:
        for contig, _, sequence in fasta_records(args.fasta):
            matches = [hap for hap, pattern in patterns.items() if pattern.search(contig)]
            if len(matches) != 1:
                excluded.write(f"{contig}\t{'unclassified' if not matches else 'matches_both_regexes'}\n")
                continue
            hap = matches[0]
            source = patterns[hap].sub("", contig).rstrip("_-.")
            write_fasta_record(fasta_handles[hap], contig, sequence)
            map_handles[hap].write(f"{source},{contig}\n")
            counts[hap] += 1
    finally:
        for handle in [*fasta_handles.values(), *map_handles.values(), excluded]:
            handle.close()
    if min(counts.values()) == 0:
        raise ValueError(
            f"could not identify both haplotypes (hap1={counts['hap1']}, hap2={counts['hap2']}); "
            "contigs must end in _hap1 or _hap2"
        )
    with open(output / "reference_summary.tsv", "w") as handle:
        handle.write("haplotype\tcontigs\n")
        for hap in ("hap1", "hap2"):
            handle.write(f"{hap}\t{counts[hap]}\n")


if __name__ == "__main__":
    main()
