#!/usr/bin/env python3
"""Gene, transcript and chromosome maps read straight from a diploid GFF3.

Only the annotation is needed: no assembly FASTA and no transcriptome are
built. This is what the ploidy route and transcript-level analyses require, and
the identifiers match the ones BUILD_ANNOTATION_BUNDLE produces.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from annotation_bundle import chromosome_label, copy_identity, first_attr
from common import parse_gff_attributes


TRANSCRIPT_TYPES = {"transcript", "mrna"}
HAPLOTYPE_RE = re.compile(r"[_-]hap([12])(?:$|[_-])", re.IGNORECASE)


def parse(path: str) -> tuple[dict[str, dict], dict[str, tuple[str, str]]]:
    """Return {gene_id: {name, hap1, hap2}} and {transcript_id: (gene_id, name)}.

    hap1 and hap2 hold (copy, chromosome) of the primary copy on each haplotype.
    """
    genes: dict[str, dict] = {}
    transcripts: dict[str, tuple[str, str]] = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 GFF3 columns")
            feature_type = fields[2].lower()
            if feature_type != "gene" and feature_type not in TRANSCRIPT_TYPES:
                continue
            attrs = parse_gff_attributes(fields[8])
            if feature_type == "gene":
                raw = first_attr(attrs, "gene_id", "ID")
                if not raw:
                    continue
                base, copy = copy_identity(raw, attrs)
                haplotype = HAPLOTYPE_RE.search(fields[0])
                gene = genes.setdefault(base, {"name": first_attr(attrs, "gene_name", "Name") or base})
                if haplotype:
                    hap = f"hap{haplotype.group(1)}"
                    gene[hap] = min(gene.get(hap, (copy, chromosome_label(fields[0]))),
                                    (copy, chromosome_label(fields[0])))
            else:
                raw_tx = first_attr(attrs, "transcript_id", "ID")
                raw_gene = first_attr(attrs, "gene_id", "Parent")
                if not raw_tx or not raw_gene:
                    continue
                base_tx, _copy = copy_identity(raw_tx, attrs)
                base_gene, _copy = copy_identity(raw_gene, attrs)
                transcripts.setdefault(base_tx, (base_gene,
                                                 first_attr(attrs, "gene_name", "Name") or base_gene))
    if not genes:
        raise ValueError(f"{path}: no gene features found")
    return genes, transcripts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    genes, transcripts = parse(args.gff3)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "gene_map.tsv", "w") as handle:
        handle.write("GID\tChrom_H1\tChrom_H2\tgene_name\n")
        for gid in sorted(genes):
            gene = genes[gid]
            chromosomes = (gene.get(hap, (0, "NA"))[1] for hap in ("hap1", "hap2"))
            handle.write(f"{gid}\t" + "\t".join(chromosomes) + f"\t{gene['name']}\n")
    with open(out / "tx2gene.tsv", "w") as handle:
        handle.write("transcript_id\tgene_id\tgene_name\n")
        for tx in sorted(transcripts):
            gid, name = transcripts[tx]
            handle.write(f"{tx}\t{gid}\t{name}\n")
    print(f"[build_feature_map] {len(genes)} genes, {len(transcripts)} transcripts")


if __name__ == "__main__":
    main()
