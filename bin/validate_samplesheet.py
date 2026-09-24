#!/usr/bin/env python3
"""Validate the single public longBayesASE-CN samplesheet format."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

from common import open_text, sniff_delimiter


COLUMNS = ["sample", "condition", "dna_id", "rna", "dna", "ploidy"]
PLOIDY_ITEM = re.compile(r"^[^:;,\s]+:(?:\d+(?:\.\d+)?):(?:\d+(?:\.\d+)?)$")


def existing_paths(value: str, base: Path, label: str, row_number: int) -> str:
    if not value:
        return ""
    resolved = []
    for item in value.split(";"):
        path = Path(os.path.expandvars(os.path.expanduser(item.strip())))
        path = path if path.is_absolute() else base / path
        path = path.resolve()
        if not path.exists():
            raise ValueError(f"row {row_number}: {label} file does not exist: {path}")
        resolved.append(str(path))
    aligned = [path.lower().endswith(".bam") for path in resolved]
    if any(aligned) and (len(resolved) != 1 or not all(aligned)):
        raise ValueError(f"row {row_number}: an aligned {label} BAM must be the only {label} file")
    return ";".join(resolved)


def validate_ploidy(value: str, row_number: int) -> str:
    value = value.strip()
    if not value:
        return ""
    items = [item.strip() for item in value.split(";") if item.strip()]
    if not items or any(not PLOIDY_ITEM.match(item) for item in items):
        raise ValueError(
            f"row {row_number}: ploidy must use chromosome:CN_H1:CN_H2;... "
            "(example: chr7:2:1;chr22:1:2)"
        )
    chromosomes = [item.split(":", 1)[0] for item in items]
    if len(chromosomes) != len(set(chromosomes)):
        raise ValueError(f"row {row_number}: ploidy contains a chromosome more than once")
    return ";".join(items)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    # diffase_priors is the diffase step when it has to simulate the priors,
    # which needs the RNA libraries.
    parser.add_argument("--step", choices=["all", "dna", "rna", "diffase", "diffase_priors"], default="all")
    args = parser.parse_args()

    source = Path(args.input).resolve()
    delimiter = sniff_delimiter(source)
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        incoming = reader.fieldnames or []
        if incoming != COLUMNS:
            raise ValueError(
                "samplesheet columns must be exactly, in this order: " + "\t".join(COLUMNS)
            )
        rows = []
        seen_samples = set()
        dna_by_id: dict[str, str] = {}
        ploidy_by_id: dict[str, str] = {}
        for row_number, raw in enumerate(reader, 2):
            row = {name: (raw.get(name) or "").strip() for name in COLUMNS}
            for required in ("sample", "condition", "dna_id"):
                if not row[required]:
                    raise ValueError(f"row {row_number}: {required} cannot be empty")
                if any(token in row[required] for token in ("/", "\\", ",", ";", "\t")):
                    raise ValueError(f"row {row_number}: invalid character in {required}: {row[required]}")
            if row["sample"] in seen_samples:
                raise ValueError(f"row {row_number}: duplicate sample: {row['sample']}")
            row["rna"] = existing_paths(row["rna"], source.parent, "RNA", row_number)
            row["dna"] = existing_paths(row["dna"], source.parent, "DNA", row_number)
            row["ploidy"] = validate_ploidy(row["ploidy"], row_number)
            if args.step in {"all", "rna"} and not row["rna"]:
                raise ValueError(f"row {row_number}: step {args.step} requires an RNA file")
            if args.step == "diffase_priors" and not row["rna"]:
                raise ValueError(f"row {row_number}: simulating the mapping priors requires an RNA file; "
                                 "give it, or pass --priors")
            if args.step == "dna" and not row["dna"]:
                raise ValueError(f"row {row_number}: step dna requires a DNA file")
            previous_dna = dna_by_id.setdefault(row["dna_id"], row["dna"])
            if row["dna"] and previous_dna and row["dna"] != previous_dna:
                raise ValueError(f"row {row_number}: dna_id {row['dna_id']} maps to different DNA files")
            previous_ploidy = ploidy_by_id.setdefault(row["dna_id"], row["ploidy"])
            if row["ploidy"] != previous_ploidy:
                raise ValueError(f"row {row_number}: dna_id {row['dna_id']} has inconsistent ploidy")
            seen_samples.add(row["sample"])
            rows.append(row)

    if not rows:
        raise ValueError("samplesheet has no data rows")
    if args.step == "all":
        dna_presence = {bool(path) for path in dna_by_id.values()}
        if len(dna_presence) > 1:
            raise ValueError("step all requires DNA for every dna_id or for none of them")
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
