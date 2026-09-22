#!/usr/bin/env python3
"""Build a tiny deterministic raw-read dataset for an end-to-end run."""

from __future__ import annotations

import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent / "toy_full"


def mutate(sequence: str, every: int, offset: int) -> str:
    bases = "ACGT"
    output = list(sequence)
    for index in range(offset, len(output), every):
        output[index] = bases[(bases.index(output[index]) + 1) % 4]
    return "".join(output)


def fastq(path: Path, libraries: dict[str, tuple[str, int]]) -> None:
    with path.open("w") as handle:
        for label, (sequence, count) in libraries.items():
            for number in range(count):
                handle.write(f"@{label}_{number + 1}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    random.seed(1701)
    base1 = "".join(random.choices("ACGT", k=800))
    gain1 = "".join(random.choices("ACGT", k=800))
    sequences = {
        "GBASE_hap1": base1,
        "GBASE_hap2": mutate(base1, 19, 3),
        "GGAIN_hap1": gain1,
        "GGAIN_hap2": mutate(gain1, 17, 5),
    }

    contigs = {
        "chr1_hap1": sequences["GBASE_hap1"],
        "chr1_hap2": sequences["GBASE_hap2"],
        "chr7_hap1": sequences["GGAIN_hap1"],
        "chr7_hap2": sequences["GGAIN_hap2"],
    }
    with (ROOT / "diploid.fa").open("w") as handle:
        for name, sequence in contigs.items():
            handle.write(f">{name}\n{sequence}\n")

    with (ROOT / "diploid.gff3").open("w") as handle:
        handle.write("##gff-version 3\n")
        for chrom, hap, gene, tx in (
            ("chr1", "hap1", "GBASE", "TBASE"),
            ("chr1", "hap2", "GBASE", "TBASE"),
            ("chr7", "hap1", "GGAIN", "TGAIN"),
            ("chr7", "hap2", "GGAIN", "TGAIN"),
        ):
            seqid = f"{chrom}_{hap}"
            attrs = f"gene_id={gene};gene_name={gene}"
            handle.write(f"{seqid}\tstub\tgene\t1\t800\t.\t+\t.\tID={gene};{attrs}\n")
            handle.write(f"{seqid}\tstub\ttranscript\t1\t800\t.\t+\t.\tID={tx};Parent={gene};{attrs};transcript_id={tx}\n")
            handle.write(f"{seqid}\tstub\texon\t1\t800\t.\t+\t.\tID={tx}.exon;Parent={tx};{attrs};transcript_id={tx}\n")

    fastq(ROOT / "DNA_WT.fastq", {name: (sequence, 12) for name, sequence in sequences.items()})
    fastq(ROOT / "DNA_A.fastq", {
        "GBASE_hap1": (sequences["GBASE_hap1"], 12),
        "GBASE_hap2": (sequences["GBASE_hap2"], 12),
        "GGAIN_hap1": (sequences["GGAIN_hap1"], 24),
        "GGAIN_hap2": (sequences["GGAIN_hap2"], 12),
    })

    rows = ["sample\tcondition\tdna_id\trna\tdna\tploidy"]
    for condition, prefix, dna_id in (("ANEUPLOID", "A", "DNA_A"), ("WT", "W", "DNA_WT")):
        for replicate in range(1, 4):
            sample = f"{prefix}{replicate}"
            h1 = 24 if condition == "ANEUPLOID" else 12
            fastq(ROOT / f"{sample}.fastq", {
                "GBASE_hap1": (sequences["GBASE_hap1"], h1),
                "GBASE_hap2": (sequences["GBASE_hap2"], 12),
                "GGAIN_hap1": (sequences["GGAIN_hap1"], h1),
                "GGAIN_hap2": (sequences["GGAIN_hap2"], 12),
            })
            gain = "chr7:2:1" if condition == "ANEUPLOID" else ""
            rows.append(
                f"{sample}\t{condition}\t{dna_id}\t{ROOT / f'{sample}.fastq'}\t"
                f"{ROOT / f'{dna_id}.fastq'}\t{gain}"
            )
    (ROOT / "samplesheet.tsv").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
