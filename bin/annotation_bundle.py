#!/usr/bin/env python3
"""Characterize two Liftoff annotations and create DNA/RNA reference tables."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from common import parse_gff_attributes


@dataclass
class GeneCopy:
    base_id: str
    feature_id: str
    gene_name: str
    hap: str
    seqid: str
    start: int
    end: int
    strand: str
    copy: int = 1
    mane: bool = False
    transcripts: set[str] = field(default_factory=set)


def copy_identity(raw_id: str, attrs: dict[str, str]) -> tuple[str, int]:
    extra = attrs.get("extra_copy_number")
    if extra not in (None, ""):
        try:
            number = int(extra) + 1
        except ValueError:
            number = 1
        base = re.sub(r"_\d+$", "", raw_id) if number > 1 else raw_id
        return base, number
    match = re.search(r"_(\d+)$", raw_id)
    if match and attrs.get("copy_num_ID"):
        return raw_id[: match.start()], int(match.group(1)) + 1
    return raw_id, 1


def first_attr(attrs: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = attrs.get(key)
        if value:
            return value.split(",")[0]
    return ""


def parse_annotation(path: str, hap: str, combined_gff):
    gene_by_raw: dict[str, GeneCopy] = {}
    genes_by_base: dict[str, list[GeneCopy]] = defaultdict(list)
    genes: list[GeneCopy] = []
    tx_to_gene_raw: dict[str, str] = {}
    tx_base: dict[str, str] = {}
    tx_stats = defaultdict(lambda: [0, 0])

    gene_by_feature: dict[str, GeneCopy] = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 GFF/GTF columns")
            combined_gff.write(line if line.endswith("\n") else line + "\n")
            attrs = parse_gff_attributes(fields[8])
            feature_type = fields[2].lower()

            if feature_type == "gene":
                raw = first_attr(attrs, "gene_id", "ID")
                if not raw:
                    continue
                base, copy = copy_identity(raw, attrs)
                gene = GeneCopy(
                    base_id=base,
                    feature_id=f"{raw}_{hap}",
                    gene_name=first_attr(attrs, "gene_name", "Name") or base,
                    hap=hap,
                    seqid=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    copy=copy,
                    mane="MANE_Select" in fields[8],
                )
                genes.append(gene)
                genes_by_base[base].append(gene)
                gene_by_feature[gene.feature_id] = gene
                for identifier in {raw, first_attr(attrs, "ID")} - {""}:
                    gene_by_raw[identifier] = gene
                continue

            if feature_type in {"transcript", "mrna"}:
                raw_tx = first_attr(attrs, "transcript_id", "ID")
                raw_gene = first_attr(attrs, "gene_id", "Parent")
                gene = gene_by_raw.get(raw_gene)
                if gene is None and raw_gene:
                    base_gene, _ = copy_identity(raw_gene, attrs)
                    candidates = genes_by_base.get(base_gene, [])
                    gene = candidates[0] if candidates else None
                if not raw_tx or gene is None:
                    continue
                base_tx, _ = copy_identity(raw_tx, attrs)
                for identifier in {raw_tx, first_attr(attrs, "ID")} - {""}:
                    tx_to_gene_raw[identifier] = gene.feature_id
                    tx_base[identifier] = base_tx
                gene.transcripts.add(base_tx)
                gene.mane = gene.mane or "MANE_Select" in fields[8]
                continue

            if feature_type == "exon":
                raw_tx = first_attr(attrs, "transcript_id", "Parent")
                feature_gene = tx_to_gene_raw.get(raw_tx)
                if feature_gene is None:
                    gene = gene_by_raw.get(first_attr(attrs, "gene_id"))
                    feature_gene = gene.feature_id if gene else None
                if feature_gene is None:
                    continue
                base_tx = tx_base.get(raw_tx) or copy_identity(raw_tx, attrs)[0]
                key = (feature_gene, base_tx, raw_tx)
                tx_stats[key][0] += 1
                # Exon length is end - start; the summary metrics are defined this way.
                tx_stats[key][1] += int(fields[4]) - int(fields[3])
                gene_by_feature[feature_gene].transcripts.add(base_tx)

    return genes, tx_stats


def chromosome_label(seqid: str) -> str:
    return re.sub(r"(?:[_-])hap[12](?:$|[_-].*$)", "", seqid, flags=re.IGNORECASE)


def primary_chromosome(copies: list[GeneCopy]) -> str:
    """Chromosome of the copy Liftoff placed first (no _N suffix), NA without copies."""
    return chromosome_label(min(copies, key=lambda gene: gene.copy).seqid) if copies else "NA"


def flag(value: bool) -> str:
    return str(value).upper()


def largest_spread(groups, index: int, relative: bool = False) -> float:
    """Largest difference between the copies of a transcript in exon number (index 0) or
    spliced length (index 1), in percent of the longest copy when relative."""
    spreads = []
    for copies in groups:
        values = [copy[index] for copy in copies]
        if values and max(values) > 0:
            spread = max(values) - min(values)
            spreads.append(100 * spread / max(values) if relative else spread)
    return max(spreads, default=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hap1-gff", required=True)
    parser.add_argument("--hap2-gff", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-length-difference", type=float, default=0.05,
                        help="largest length difference between the copies of a transcript, as a fraction "
                             "of its length, for a single-copy gene to be a Singleton")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "diploid.gff3", "w") as combined_gff:
        combined_gff.write("##gff-version 3\n")
        parsed = {
            "hap1": parse_annotation(args.hap1_gff, "hap1", combined_gff),
            "hap2": parse_annotation(args.hap2_gff, "hap2", combined_gff),
        }

    all_genes = [gene for hap in parsed.values() for gene in hap[0]]
    genes_by_base_hap = defaultdict(list)
    for gene in all_genes:
        genes_by_base_hap[(gene.base_id, gene.hap)].append(gene)
    all_ids = sorted({gene.base_id for gene in all_genes})

    tx_by_gene_hap = defaultdict(lambda: defaultdict(list))
    for hap, (hap_genes, stats) in parsed.items():
        feature_to_base = {gene.feature_id: gene.base_id for gene in hap_genes}
        for (feature_gene, base_tx, _raw_tx), (n_exons, length) in stats.items():
            tx_by_gene_hap[(feature_to_base[feature_gene], hap)][base_tx].append((n_exons, length))

    summary_rows = []
    for gid in all_ids:
        h1 = genes_by_base_hap.get((gid, "hap1"), [])
        h2 = genes_by_base_hap.get((gid, "hap2"), [])
        tx_h1 = tx_by_gene_hap.get((gid, "hap1"), {})
        tx_h2 = tx_by_gene_hap.get((gid, "hap2"), {})
        all_tx = set(tx_h1) | set(tx_h2)
        both = [tx_h1.get(tx, []) + tx_h2.get(tx, []) for tx in all_tx]
        exons_diff = largest_spread(both, 0)
        length_diff_pct = largest_spread(both, 1, relative=True)
        every_tx_both = bool(all_tx) and all(tx in tx_h1 and tx in tx_h2 for tx in all_tx)
        # Fewest haplotypes any transcript of this haplotype is annotated on (0 without transcripts).
        min_tx = {hap: min((1 + (tx in other) for tx in txs), default=0)
                  for hap, txs, other in (("hap1", tx_h1, tx_h2), ("hap2", tx_h2, tx_h1))}
        # Copies of the gene on each chromosome of each haplotype.
        layout_h1 = Counter(chromosome_label(g.seqid) for g in h1)
        layout_h2 = Counter(chromosome_label(g.seqid) for g in h2)
        single_copy = len(h1) == len(h2) == 1
        structural = not every_tx_both or exons_diff > 0 or length_diff_pct > 100 * args.max_length_difference
        if not (h1 and h2):
            gene_class = "One_haplotype"
        elif single_copy and layout_h1 != layout_h2:
            gene_class = "Single_copy_different_chromosomes"
        elif single_copy:
            gene_class = "Structural_difference" if structural else "Singleton"
        elif layout_h1 == layout_h2:
            gene_class = "Paired_copies"
        else:
            several = len(layout_h1 | layout_h2) > 1
            gene_class = "Unpaired_copies_several_chromosomes" if several else "Unpaired_copies_same_chromosome"
        gene_name = (h1 or h2)[0].gene_name
        summary_rows.append(
            {
                "GID": gid,
                "gene_name": gene_name,
                "Chrom_H1": primary_chromosome(h1),
                "Chrom_H2": primary_chromosome(h2),
                "copies_h1": len(h1),
                "copies_h2": len(h2),
                "chromosomes_h1": ",".join(sorted(layout_h1)) or "NA",
                "chromosomes_h2": ",".join(sorted(layout_h2)) or "NA",
                "length_diff_h1": largest_spread(tx_h1.values(), 1),
                "length_diff_h2": largest_spread(tx_h2.values(), 1),
                "length_diff_haps": largest_spread(both, 1),
                "length_diff_pct_haps": round(length_diff_pct, 2),
                "exons_diff_h1": largest_spread(tx_h1.values(), 0),
                "exons_diff_h2": largest_spread(tx_h2.values(), 0),
                "exons_diff_haps": exons_diff,
                "min_tx_copies_h1": min_tx["hap1"],
                "min_tx_copies_h2": min_tx["hap2"],
                "MANE_h1": flag(any(g.mane for g in h1)),
                "MANE_h2": flag(any(g.mane for g in h2)),
                "MANE_tag": "MANE" if any(g.mane for g in h1 + h2) else "Non_MANE",
                "Class": gene_class,
            }
        )

    summary_fields = list(summary_rows[0]) if summary_rows else ["GID"]
    with open(out / "annotation_summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    sorted_genes = sorted(all_genes, key=lambda g: (g.seqid, g.start, g.end, g.feature_id))
    with open(out / "genes.bed", "w") as bed, open(out / "genes.saf", "w") as saf:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for gene in sorted_genes:
            bed.write(f"{gene.seqid}\t{max(0, gene.start - 1)}\t{gene.end}\t{gene.feature_id}\n")
            saf.write(f"{gene.feature_id}\t{gene.seqid}\t{gene.start}\t{gene.end}\t{gene.strand}\n")

    tx_gene = {}
    for hap, (genes, stats) in parsed.items():
        feature_to_base = {gene.feature_id: gene.base_id for gene in genes}
        feature_to_name = {gene.feature_id: gene.gene_name for gene in genes}
        for (feature_gene, base_tx, _), _values in stats.items():
            tx_gene.setdefault(base_tx, (feature_to_base[feature_gene], feature_to_name[feature_gene]))
    with open(out / "tx2gene.tsv", "w") as handle:
        handle.write("transcript_id\tgene_id\tgene_name\n")
        for tx in sorted(tx_gene):
            gid, name = tx_gene[tx]
            handle.write(f"{tx}\t{gid}\t{name}\n")


if __name__ == "__main__":
    main()
