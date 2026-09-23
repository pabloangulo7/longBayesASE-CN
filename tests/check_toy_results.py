#!/usr/bin/env python3
"""Assert the toy dataset behaves as designed. Used by CI, runnable by hand."""

import csv
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "test-results/diffase/diffASE_results.tsv"
rows = {row["ID"]: row for row in csv.DictReader(open(path), delimiter="\t")}

assert set(rows) == {"GENE_BALANCED", "GENE_DIFFASE"}, f"unexpected genes: {sorted(rows)}"
valid_fits = {"Success", "pvalue0"}
assert all(row["analysis_flag"] in valid_fits for row in rows.values()), rows

# GENE_BALANCED doubles H1 exactly with its copy number, so once the copy-number term is
# in the model it is not differential. GENE_DIFFASE keeps the duplicated
# haplotype flat, which is allele-specific dosage compensation, so it is.
balanced = float(rows["GENE_BALANCED"]["diffAI_pvalue"])
differential = float(rows["GENE_DIFFASE"]["diffAI_pvalue"])
assert differential < 0.05, f"GENE_DIFFASE should be differential, got p={differential}"
assert balanced > 0.05, f"GENE_BALANCED scales with copy number, got p={balanced}"
print(f"OK: balanced p={balanced:.3f}, differential p={differential:.3f}")
