#!/usr/bin/env python3
"""Keep one transcript copy per haplotype, label it, and concatenate haplotypes."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from common import fasta_records, write_fasta_record


def load_haplotype(path: str, hap: str) -> tuple[dict[str, str], list[str]]:
    selected: dict[str, str] = {}
    order: list[str] = []
    for raw_id, _header, sequence in fasta_records(path):
        # Extra Liftoff copies end in _<number>; only the primary copy is kept.
        if re.search(r"_\d+$", raw_id):
            continue
        if raw_id not in selected:
            order.append(raw_id)
            selected[raw_id] = sequence
    return {f"{base}_{hap}": selected[base] for base in order}, order


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    datasets = {}
    for hap, path in (("hap1", args.hap1), ("hap2", args.hap2)):
        records, order = load_haplotype(path, hap)
        datasets[hap] = records
        with open(output / f"{hap}_transcriptome.fa", "w") as handle:
            for base in order:
                identifier = f"{base}_{hap}"
                write_fasta_record(handle, identifier, records[identifier])
    with open(output / "diploid_transcriptome.fa", "w") as handle:
        for hap in ("hap1", "hap2"):
            for identifier, sequence in datasets[hap].items():
                write_fasta_record(handle, identifier, sequence)
    with open(output / "transcriptome_summary.tsv", "w") as handle:
        handle.write("haplotype\ttranscripts\n")
        for hap in ("hap1", "hap2"):
            handle.write(f"{hap}\t{len(datasets[hap])}\n")


if __name__ == "__main__":
    main()
