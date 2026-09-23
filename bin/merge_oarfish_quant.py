#!/usr/bin/env python3
"""Collapse Oarfish estimates over haplotypes into gene and transcript tables.

Oarfish quantifies every transcript of the diploid transcriptome, so each
transcript appears once per haplotype. How the EM splits a read between the two
haplotypes of a transcript depends on the haplotype abundances it estimated, so
those per-haplotype values are not used; their sum is the transcript's
expression, and the sum over a gene's transcripts is the gene's. Both tables are
long, with the columns sample, ID and num_reads.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict

from common import files_by_sample, gene_of, load_tx2gene, split_haplotype


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quant", nargs="+", required=True, help="per-sample SAMPLE.quant files")
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--genes-output", required=True)
    parser.add_argument("--transcripts-output", required=True)
    args = parser.parse_args()

    tx2gene = load_tx2gene(args.tx2gene)
    by_sample = files_by_sample(args.quant, r"\.quant")

    with open(args.genes_output, "w") as genes_out, open(args.transcripts_output, "w") as tx_out:
        genes_out.write("sample\tID\tnum_reads\n")
        tx_out.write("sample\tID\tnum_reads\n")
        for sample in sorted(by_sample):
            transcripts: dict[str, float] = defaultdict(float)
            with open(by_sample[sample]) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    transcript, _hap = split_haplotype(row["tname"])
                    transcripts[transcript] += float(row["num_reads"])
            genes: dict[str, float] = defaultdict(float)
            for transcript, value in transcripts.items():
                genes[gene_of(tx2gene, transcript)] += value
            for transcript in sorted(transcripts):
                tx_out.write(f"{sample}\t{transcript}\t{transcripts[transcript]:.3f}\n")
            for gene in sorted(genes):
                genes_out.write(f"{sample}\t{gene}\t{genes[gene]:.3f}\n")


if __name__ == "__main__":
    main()
