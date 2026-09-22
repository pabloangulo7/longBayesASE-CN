#!/usr/bin/env python3
"""Convert Oarfish read-to-transcript probabilities into a simple TSV."""

from __future__ import annotations

import argparse
import csv
import re


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prob", required=True)
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.tx2gene) as handle:
        tx_map = {row["transcript_id"]: row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
    with open(args.prob) as source, open(args.output, "w") as target:
        header = source.readline().strip().split()
        if not header:
            raise ValueError("empty Oarfish probability file")
        transcript_count = int(header[0])
        transcripts = [source.readline().strip().split()[0] for _ in range(transcript_count)]
        target.write("read_id\ttranscript_id\tgene_id\tprobability\n")
        for line_number, line in enumerate(source, transcript_count + 2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                fields = line.split()
            if len(fields) < 2:
                continue
            read_id = fields[0]
            assignments = int(fields[1])
            if len(fields) < 2 + 2 * assignments:
                raise ValueError(f"malformed Oarfish probability row at line {line_number}")
            indices = [int(value) for value in fields[2 : 2 + assignments]]
            probabilities = [float(value) for value in fields[2 + assignments : 2 + 2 * assignments]]
            for index, probability in zip(indices, probabilities):
                transcript = transcripts[index]
                match = re.search(r"_hap([12])$", transcript)
                if not match:
                    raise ValueError(f"transcript lacks _hap1/_hap2 suffix: {transcript}")
                hap = f"hap{match.group(1)}"
                base_tx = transcript[: match.start()]
                lookup_tx = base_tx.split("|", 1)[0]
                gene = tx_map.get(base_tx) or tx_map.get(lookup_tx)
                if gene is None:
                    raise ValueError(f"no tx2gene entry for transcript: {base_tx}")
                target.write(f"{read_id}\t{transcript}\t{gene}_{hap}\t{probability:.12g}\n")


if __name__ == "__main__":
    main()
