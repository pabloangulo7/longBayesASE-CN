#!/usr/bin/env python3
"""Shared parsing helpers for longBayesASE-CN command-line utilities."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, TextIO, Tuple
from urllib.parse import unquote


HAP_SUFFIX_RE = re.compile(r"_hap([12])$")
COPY_SUFFIX_RE = re.compile(r"_\d+$")

# A haplotype copy number of exactly 0 makes the negative binomial mean 0 and
# the sampler fails, so a lost haplotype is floored here.
MIN_COPY_NUMBER = 0.05


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def sniff_delimiter(path: str | Path) -> str:
    with open_text(path) as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
    except csv.Error:
        return "\t" if "\t" in sample else ","


def parse_gff_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for field in raw.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def feature_parts(feature_id: str) -> Tuple[str, str | None]:
    """Return (base feature, H1/H2) from a pipeline feature ID such as GENE_1_hap2."""
    hap_match = HAP_SUFFIX_RE.search(feature_id)
    hap = f"H{hap_match.group(1)}" if hap_match else None
    base = COPY_SUFFIX_RE.sub("", HAP_SUFFIX_RE.sub("", feature_id))
    return base, hap


def split_haplotype(name: str) -> Tuple[str, int]:
    """Return (transcript, 1 or 2) from a diploid transcriptome name such as TX1_hap2."""
    match = HAP_SUFFIX_RE.search(name)
    if not match:
        raise ValueError(f"transcript lacks _hap1/_hap2 suffix: {name}")
    return name[: match.start()], int(match.group(1))


def gene_of(tx2gene: Dict[str, str], transcript: str) -> str:
    gene = tx2gene.get(transcript)
    if gene is None:
        raise ValueError(f"no tx2gene entry for transcript: {transcript}")
    return gene


def files_by_sample(paths: Iterable[str], suffix: str) -> Dict[str, str]:
    """Map each sample to its per-library file, named <sample><suffix>.

    suffix is a regular expression; two files for one sample are an error.
    """
    pattern = re.compile(rf"^(?P<sample>.+){suffix}$")
    by_sample: Dict[str, str] = {}
    for path in paths:
        match = pattern.match(Path(path).name)
        if not match:
            raise ValueError(f"cannot identify the sample from the file name: {path}")
        sample = match.group("sample")
        if sample in by_sample:
            raise ValueError(f"two files for sample {sample}: {by_sample[sample]} and {path}")
        by_sample[sample] = path
    return by_sample


def fasta_records(path: str | Path) -> Iterator[Tuple[str, str, str]]:
    """Yield FASTA identifier, complete header, and sequence."""
    ident = ""
    header = ""
    sequence: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    yield ident, header, "".join(sequence)
                header = line[1:].rstrip("\n")
                ident = header.split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if header:
        yield ident, header, "".join(sequence)


def write_fasta_record(handle: TextIO, identifier: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def load_tx2gene(path: str | Path | None) -> Dict[str, str]:
    """Map transcript_id to gene_id; an absent or header-only table gives {}."""
    if not path:
        return {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"transcript_id", "gene_id"}.issubset(reader.fieldnames or []):
            return {}
        return {row["transcript_id"]: row["gene_id"] for row in reader}


def parse_ploidy(value: str) -> Dict[str, Tuple[float, float]]:
    """Parse the samplesheet ploidy column, chromosome:CN_H1:CN_H2;..."""
    result = {}
    for item in (value or "").split(";"):
        if item:
            chromosome, h1, h2 = item.split(":")
            result[chromosome] = float(h1), float(h2)
    return result


def median(values: Iterable[float]) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return float("nan")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0
