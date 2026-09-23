#!/usr/bin/env python3
"""Check the raw-read two-gene end-to-end test."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "full-integration-results"


def read(path: Path):
    with path.open() as handle:
        rows = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(rows, delimiter="\t"))


qc = read(ROOT / "dna" / "copy_number.qc.tsv")
gain = next(row for row in qc if row["sample"] == "DNA_A" and row["Chrom"] == "chr7")
assert gain["changed"] == "True"
assert gain["direction"] == "gain"
assert abs(float(gain["CN_H1"]) - 2) < 0.05
assert abs(float(gain["CN_H2"]) - 1) < 0.05

priors = read(ROOT / "priors" / "mapping_priors.gene.tsv")
assert {row["ID"] for row in priors} == {"GBASE", "GGAIN"}

results = {row["ID"]: row for row in read(ROOT / "diffase" / "diffASE_results.tsv")}
required = {
    "groupA_alphaAI_pvalue", "groupA_thetaAI_pvalue",
    "groupB_alphaAI_pvalue", "groupB_thetaAI_pvalue",
    "diffAI_pvalue", "rope_value", "analysis_flag",
}
assert required.issubset(results["GBASE"])
assert results["GBASE"]["analysis_flag"] in {"Success", "pvalue0"}
assert results["GGAIN"]["analysis_flag"] in {"Success", "pvalue0"}
assert float(results["GBASE"]["diffAI_pvalue"]) < 0.1
assert float(results["GGAIN"]["diffAI_pvalue"]) > 0.1
assert float(results["GGAIN"]["rope_value"]) > 0.5
