#!/usr/bin/env python3
"""Split a diploid GFF3 by the required _hap1/_hap2 contig suffixes."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    args = parser.parse_args()

    counts = {"hap1": 0, "hap2": 0}
    with open(args.gff3) as source, open(args.hap1, "w") as h1, open(args.hap2, "w") as h2:
        h1.write("##gff-version 3\n")
        h2.write("##gff-version 3\n")
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{args.gff3}:{line_number}: expected 9 GFF3 columns")
            if "_hap1" in fields[0]:
                h1.write(line if line.endswith("\n") else line + "\n")
                counts["hap1"] += 1
            elif "_hap2" in fields[0]:
                h2.write(line if line.endswith("\n") else line + "\n")
                counts["hap2"] += 1
            else:
                raise ValueError(
                    f"{args.gff3}:{line_number}: contig {fields[0]} does not contain _hap1 or _hap2"
                )
    if min(counts.values()) == 0:
        raise ValueError("the GFF3 must contain features on both _hap1 and _hap2 contigs")


if __name__ == "__main__":
    main()
