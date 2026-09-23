from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pysam
import pytest


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"


def run_script(name: str, *args: object) -> None:
    subprocess.run([sys.executable, str(BIN / name), *map(str, args)], check=True)


def read_tsv(path: Path):
    with path.open() as handle:
        rows = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(rows, delimiter="\t"))


TRANSCRIPTS = ["T1_hap1", "T1_hap2", "T1B_hap1", "T1B_hap2", "T2_hap1", "T2_hap2"]


def write_rna_bam(path: Path, alignments: list[tuple[str, str, int, int]]) -> None:
    """Write (read, transcript, AS, flag) records in the given, unsorted order."""
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": name, "LN": 1000} for name in TRANSCRIPTS]}
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for read, transcript, score, flag in alignments:
            record = pysam.AlignedSegment()
            record.query_name = read
            record.flag = flag
            record.reference_id = TRANSCRIPTS.index(transcript)
            record.reference_start = 0
            record.mapping_quality = 60
            record.cigar = ((0, 50),)
            if not flag & 256:
                record.query_sequence = "A" * 50
                record.query_qualities = pysam.qualitystring_to_array("I" * 50)
            record.set_tag("AS", score)
            handle.write(record)


def rna_tx2gene(tmp_path: Path) -> Path:
    path = tmp_path / "tx2gene.tsv"
    path.write_text("transcript_id\tgene_id\tgene_name\nT1\tG1\tG1\nT1B\tG1\tG1\nT2\tG2\tG2\n")
    return path


def test_rna_haplotypes_come_from_the_best_alignments(tmp_path: Path) -> None:
    bam = tmp_path / "reads.bam"
    write_rna_bam(bam, [
        ("only_hap1", "T1_hap1", 100, 0),
        ("both_haps", "T1_hap1", 100, 0), ("both_haps", "T1_hap2", 100, 256),
        # A lower score on hap2 is sequence evidence for hap1.
        ("snp", "T1_hap1", 100, 0), ("snp", "T1_hap2", 90, 256),
        ("two_isoforms", "T1_hap1", 100, 0), ("two_isoforms", "T1B_hap1", 100, 256),
        ("two_genes", "T1_hap1", 100, 0), ("two_genes", "T2_hap2", 100, 256),
        ("antisense", "T1_hap2", 100, 16),
        ("chimera", "T1_hap1", 100, 0), ("chimera", "T2_hap1", 100, 2048),
    ])

    def run(strand: str) -> dict[str, dict[str, str]]:
        prefix = tmp_path / strand
        run_script("count_rna_haplotypes.py", "--bam", bam, "--tx2gene", rna_tx2gene(tmp_path),
                   "--strand", strand, "--output-prefix", prefix)
        return {level: {row["Read_ID"]: row["Group"] for row in read_tsv(Path(f"{prefix}.{level}.readgroups.tsv"))}
                for level in ("genes", "isoforms")}

    groups = run("fw")
    assert groups["genes"] == {
        "only_hap1": "H1", "both_haps": "NonHS", "snp": "H1", "two_isoforms": "H1_multimapping",
        "two_genes": "NonHS_complex", "chimera": "H1",
    }
    assert groups["isoforms"]["two_isoforms"] == "H1_multimapping_multigene"
    qc = {row["metric"]: row["reads"] for row in read_tsv(tmp_path / "fw.genes.qc.tsv")}
    assert qc["Discarded_orientation"] == "1"
    counts = {row["ID"]: row for row in read_tsv(tmp_path / "fw.genes.counts.tsv")}
    assert counts["G1"]["H1"] == "3" and counts["G1"]["NonHS"] == "1"
    assert run("both")["genes"]["antisense"] == "H2"


def test_rna_counting_rejects_a_read_split_into_two_blocks(tmp_path: Path) -> None:
    """Counting streams one read at a time, so a read must not reappear later."""
    bam = tmp_path / "sorted.bam"
    write_rna_bam(bam, [("r1", "T1_hap1", 100, 0), ("r2", "T1_hap2", 100, 0), ("r1", "T1_hap2", 100, 256)])
    with pytest.raises(subprocess.CalledProcessError):
        run_script("count_rna_haplotypes.py", "--bam", bam, "--tx2gene", rna_tx2gene(tmp_path),
                   "--output-prefix", tmp_path / "split")


def test_oarfish_estimates_are_summed_over_haplotypes(tmp_path: Path) -> None:
    quant = []
    for sample, values in (("S1", (10.5, 4.5, 3, 0, 7, 1)), ("S2", (1, 1, 0, 0, 0, 0))):
        path = tmp_path / f"{sample}.quant"
        path.write_text("tname\tlen\tnum_reads\n" + "".join(
            f"{name}\t1000\t{value}\n" for name, value in zip(TRANSCRIPTS, values)))
        quant.append(path)
    run_script("merge_oarfish_quant.py", "--quant", *quant, "--tx2gene", rna_tx2gene(tmp_path),
               "--genes-output", tmp_path / "genes.tsv", "--transcripts-output", tmp_path / "tx.tsv")
    genes = {(row["sample"], row["ID"]): float(row["num_reads"]) for row in read_tsv(tmp_path / "genes.tsv")}
    transcripts = {(row["sample"], row["ID"]): float(row["num_reads"]) for row in read_tsv(tmp_path / "tx.tsv")}
    assert transcripts[("S1", "T1")] == 15 and transcripts[("S1", "T1B")] == 3
    assert genes == {("S1", "G1"): 18, ("S1", "G2"): 8, ("S2", "G1"): 2, ("S2", "G2"): 0}


def test_mapping_priors(tmp_path: Path) -> None:
    header = "ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\n"
    h1 = tmp_path / "h1.tsv"
    h2 = tmp_path / "h2.tsv"
    h1.write_text(header + "G1\t80\t0\t0\t0\t20\t0\n")
    h2.write_text(header + "G1\t0\t0\t75\t0\t25\t0\n")
    output = tmp_path / "priors.tsv"
    run_script("build_priors.py", "--hap1-counts", h1, "--hap2-counts", h2, "--output", output)
    row = read_tsv(output)[0]
    assert float(row["H1_prior"]) == 0.8
    assert float(row["H2_prior"]) == 0.75


def test_dna_copy_multimapping_classification(tmp_path: Path) -> None:
    featurecounts = tmp_path / "featurecounts.tsv"
    featurecounts.write_text(
        "r1\tAssigned\tchr1\tG1_1_hap1\n"
        "r1\tAssigned\tchr1\tG1_2_hap1\n"
    )
    readgroups = tmp_path / "readgroups.tsv"
    readgroups.write_text("Read_ID\tGroup\nr1\tH1_multimapping\n")
    prefix = tmp_path / "dna"
    run_script(
        "classify_dna_genes.py",
        "--featurecounts", featurecounts,
        "--readgroups", readgroups,
        "--output-prefix", prefix,
    )
    groups = read_tsv(tmp_path / "dna.readgroups_genes.tsv")
    assert groups == [{"Read_ID": "r1", "Group": "H1_multimapping_multigene_copy"}]
    counts = read_tsv(tmp_path / "dna.gene_counts.tsv")
    assert counts[0]["H1_multimapping_multigene_copy"] == "1"


def test_dna_split_uses_coordinate_bam_and_readgroup_dictionary(tmp_path: Path) -> None:
    bam = tmp_path / "maxscore.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1_hap1", "LN": 1000},
            {"SN": "chr1_hap2", "LN": 1000},
            {"SN": "chr2_hap1", "LN": 1000},
        ],
    }

    def alignment(read: str, reference: int, start: int) -> pysam.AlignedSegment:
        record = pysam.AlignedSegment()
        record.query_name = read
        record.query_sequence = "A" * 50
        record.flag = 0
        record.reference_id = reference
        record.reference_start = start
        record.mapping_quality = 60
        record.cigar = ((0, 50),)
        record.query_qualities = pysam.qualitystring_to_array("I" * 50)
        return record

    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        # r1 has two co-optimal H1 loci: the first one in coordinate order,
        # chr1_hap1, is kept.
        handle.write(alignment("r1", 0, 100))
        handle.write(alignment("r2", 0, 200))
        handle.write(alignment("r2", 1, 200))
        handle.write(alignment("r1", 2, 10))

    groups = tmp_path / "groups.tsv"
    groups.write_text("Read_ID\tGroup\nr2\tNonHS\nr1\tH1_multimapping\n")
    prefix = tmp_path / "split"
    run_script(
        "split_dna_bam.py",
        "--bam", bam,
        "--readgroups", groups,
        "--output-prefix", prefix,
    )
    with pysam.AlignmentFile(tmp_path / "split.HS.bam", "rb") as handle:
        kept = list(handle)
        assert len(kept) == 1
        assert kept[0].reference_name == "chr1_hap1"
    with pysam.AlignmentFile(tmp_path / "split.NonHS.bam", "rb") as handle:
        assert [record.reference_name for record in handle] == ["chr1_hap1", "chr1_hap2"]


def test_balanced_ase_shards(tmp_path: Path) -> None:
    output = tmp_path / "shards"
    run_script(
        "split_ase_input.py",
        "--input", ROOT / "tests/data/ase_input.tsv",
        "--contrast", "ANEUPLOID:WT",
        "--shards", 2,
        "--output-dir", output,
    )
    manifest = read_tsv(output / "manifest.tsv")
    assert len(manifest) == 2
    assert sum(int(row["genes"]) for row in manifest) == 2
    assert all(Path(output / row["file"]).exists() for row in manifest)
    shard_rows = read_tsv(output / manifest[0]["file"])
    assert {"group", "test", "groupA", "groupB"}.issubset(shard_rows[0])


def test_reference_split_uses_required_haplotype_suffixes(tmp_path: Path) -> None:
    fasta = tmp_path / "diploid.fa"
    fasta.write_text(">chr1_hap1\nACGT\n>chr1_hap2\nTGCA\n")
    out = tmp_path / "reference"
    run_script("prepare_reference.py", "--fasta", fasta, "--output-dir", out)
    assert (out / "hap1.fa").read_text().startswith(">chr1_hap1")
    assert (out / "hap2.fa").read_text().startswith(">chr1_hap2")
    assert (out / "hap1.chroms.csv").read_text() == "chr1,chr1_hap1\n"


def test_existing_diploid_gff_is_split_by_contig(tmp_path: Path) -> None:
    gff = tmp_path / "diploid.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1_hap1\ttest\tgene\t1\t10\t.\t+\t.\tID=G1;gene_id=G1\n"
        "chr1_hap2\ttest\tgene\t1\t10\t.\t+\t.\tID=G1;gene_id=G1\n"
    )
    hap1, hap2 = tmp_path / "hap1.gff3", tmp_path / "hap2.gff3"
    run_script("split_diploid_gff.py", "--gff3", gff, "--hap1", hap1, "--hap2", hap2)
    assert "chr1_hap1" in hap1.read_text()
    assert "chr1_hap2" not in hap1.read_text()
    assert "chr1_hap2" in hap2.read_text()


def test_simple_samplesheet_and_ploidy(tmp_path: Path) -> None:
    source = tmp_path / "samplesheet.tsv"
    source.write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\n"
        "Sample2_R1\tTreatment\tSample2\t\t\tchr7:2:1;chr22:1:2\n"
    )
    output = tmp_path / "validated.tsv"
    run_script(
        "validate_samplesheet.py",
        "--input", source,
        "--output", output,
        "--step", "diffase",
    )
    row = read_tsv(output)[0]
    assert row["condition"] == "Treatment"
    assert row["dna_id"] == "Sample2"
    assert row["ploidy"] == "chr7:2:1;chr22:1:2"


def test_samplesheet_accepts_multiple_technical_read_files(tmp_path: Path) -> None:
    dna1 = tmp_path / "dna_rep1.fastq"
    dna2 = tmp_path / "dna_rep2.fastq"
    rna = tmp_path / "rna.fastq"
    for path in (dna1, dna2, rna):
        path.write_text("@r\nACGT\n+\nIIII\n")
    source = tmp_path / "samples.tsv"
    source.write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\n"
        f"S1\tWT\tD1\t{rna.name}\t{dna1.name};{dna2.name}\t\n"
    )
    output = tmp_path / "validated.tsv"
    run_script("validate_samplesheet.py", "--input", source, "--output", output, "--step", "all")
    row = read_tsv(output)[0]
    assert row["dna"].split(";") == [str(dna1.resolve()), str(dna2.resolve())]


def test_unified_rna_counts(tmp_path: Path) -> None:
    counts = tmp_path / "S1.genes.counts.tsv"
    counts.write_text(
        "ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\n"
        "G1\t10\t2\t8\t1\t4\t3\n"
    )
    output = tmp_path / "gene_counts.tsv"
    run_script("merge_rna_counts.py", "--counts", counts, "--output", output)
    assert read_tsv(output) == [{"sample": "S1", "ID": "G1", "H1": "12", "H2": "9", "NonHS": "7"}]


def test_copy_number_from_samplesheet_ploidy(tmp_path: Path) -> None:
    samplesheet = tmp_path / "samples.tsv"
    samplesheet.write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\n"
        "S1\tA\tD1\t\t\tchr7:2:1;chr9:1:0\n"
    )
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text("GID\tChrom\nG1\tchr1\nG7\tchr7\nG9\tchr9\n")
    output = tmp_path / "copy.tsv"
    run_script(
        "copy_number_from_ploidy.py", "--samplesheet", samplesheet,
        "--annotation-summary", annotation, "--output", output,
    )
    rows = {row["ID"]: row for row in read_tsv(output)}
    assert (rows["G1"]["CN_H1"], rows["G1"]["CN_H2"]) == ("1.0", "1.0")
    assert (rows["G7"]["CN_H1"], rows["G7"]["CN_H2"]) == ("2.0", "1.0")
    # A declared loss is floored like a measured one: a copy number of 0 breaks the model.
    assert (rows["G9"]["CN_H1"], rows["G9"]["CN_H2"]) == ("1.0", "0.05")


def test_exhaustive_transcript_tiling(tmp_path: Path) -> None:
    fasta = tmp_path / "tx.fa"
    fasta.write_text(">TX1_hap1\nAACCGGTTAA\n>TX2_hap1\nACGT\n")
    output = tmp_path / "tiles.fastq.gz"
    run_script(
        "tile_transcripts.py",
        "--transcriptome", fasta,
        "--read-length", 6,
        "--step", 4,
        "--output", output,
    )
    with gzip.open(output, "rt") as handle:
        lines = handle.read().splitlines()
    assert lines[0::4] == [
        "@TX1_hap1__tile_1_6",
        "@TX1_hap1__tile_5_10",
        "@TX2_hap1__tile_1_4",
    ]
    assert lines[1::4] == ["AACCGG", "GGTTAA", "ACGT"]


def test_chromosome_copy_number_and_haplotype_proportions(tmp_path: Path) -> None:
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text(
        "GID\tChrom\tSet\n"
        "G0\tchr1\tSingleton\n"
        "G1\tchr1\tSingleton\n"
        "G7A\tchr7\tSingleton\n"
        "G7B\tchr7\tSingleton\n"
        "GCOPY\tchr7\tMulticopy\n"
    )
    samplesheet = tmp_path / "samplesheet.tsv"
    samplesheet.write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\nS1\tA\tD1\t\t\tchr7:2:1\n"
    )
    coverage_rows = {
        "HS": (
            "G0_hap1\t10\nG0_hap2\t10\nG1_hap1\t10\nG1_hap2\t10\n"
            "G7A_hap1\t20\nG7A_hap2\t10\nG7B_hap1\t20\nG7B_hap2\t10\n"
            "GCOPY_hap1\t3\nGCOPY_hap2\t1\n"
        ),
        "NonHS": (
            "G0_hap1\t0\nG0_hap2\t0\nG1_hap1\t0\nG1_hap2\t0\n"
            "G7A_hap1\t0\nG7A_hap2\t0\nG7B_hap1\t0\nG7B_hap2\t0\n"
            "GCOPY_hap1\t0\nGCOPY_hap2\t0\n"
        ),
        "HS_copy": "GCOPY_1_hap1\t6\nGCOPY_1_hap2\t0\n",
        "NonHS_copy": "GCOPY_1_hap1\t0\nGCOPY_1_hap2\t0\n",
    }
    coverage = []
    for group, content in coverage_rows.items():
        path = tmp_path / f"D1.{group}.coverage.tsv"
        path.write_text(content)
        coverage.append(path)
    output = tmp_path / "copy_number.tsv"
    qc = tmp_path / "copy_number.qc.tsv"
    run_script(
        "compute_copy_number.py",
        "--samplesheet", samplesheet,
        "--annotation-summary", annotation,
        "--coverage", *coverage,
        "--output", output,
        "--qc-output", qc,
    )
    rows = {row["ID"]: row for row in read_tsv(output)}
    assert float(rows["G0"]["CN_H1"]) == 1.0
    assert float(rows["G0"]["CN_H2"]) == 1.0
    assert float(rows["G7A"]["CN_H1"]) == 2.0
    assert float(rows["G7A"]["CN_H2"]) == 1.0
    assert rows["GCOPY"]["CN_H1"] == rows["G7A"]["CN_H1"]
    assert rows["GCOPY"]["CN_H2"] == rows["G7A"]["CN_H2"]
    qc_rows = {row["Chrom"]: row for row in read_tsv(qc)}
    assert float(qc_rows["chr1"]["Chrom_CN"]) == 2.0
    assert float(qc_rows["chr7"]["Chrom_CN"]) == 3.0
    assert float(qc_rows["chr7"]["H1_prop"]) == 2 / 3

    gene_output = tmp_path / "copy_number.gene.tsv"
    run_script(
        "compute_copy_number.py",
        "--samplesheet", samplesheet,
        "--annotation-summary", annotation,
        "--coverage", *coverage,
        "--cn-resolution", "gene",
        "--output", gene_output,
        "--qc-output", tmp_path / "copy_number.gene.qc.tsv",
    )
    gene_rows = {row["ID"]: row for row in read_tsv(gene_output)}
    assert abs(float(gene_rows["GCOPY"]["CN_H1"]) - 0.9) < 1e-12
    assert abs(float(gene_rows["GCOPY"]["CN_H2"]) - 0.1) < 1e-12


def test_tiling_unweighted_is_unchanged(tmp_path: Path) -> None:
    fasta = tmp_path / "tx.fa"
    fasta.write_text(">TX1_hap1\n" + "ACGT" * 400 + "\n>TX2_hap1\n" + "ACGT" * 100 + "\n")
    plain = tmp_path / "plain.fastq"
    run_script("tile_transcripts.py", "--transcriptome", fasta, "--read-length", "500",
               "--step", "250", "--output", plain)
    empty = tmp_path / "usage.tsv"
    empty.write_text("sample\tID\tnum_reads\n")
    with_flag = tmp_path / "flagged.fastq"
    run_script("tile_transcripts.py", "--transcriptome", fasta, "--read-length", "500",
               "--step", "250", "--isoform-usage", empty, "--max-replicates", "3",
               "--output", with_flag)
    assert plain.read_text() == with_flag.read_text()


def test_tiling_weights_by_isoform_usage(tmp_path: Path) -> None:
    """Tiles are replicated by the Oarfish estimates summed over every sample."""
    fasta = tmp_path / "tx.fa"
    fasta.write_text(">TX1_hap1\n" + "ACGT" * 400 + "\n>TX2_hap1\n" + "ACGT" * 400 + "\n")
    usage = tmp_path / "transcript_quant.tsv"
    usage.write_text("sample\tID\tnum_reads\nS1\tTX1\t99.5\nS1\tTX2\t0.5\nS2\tTX2\t0.5\n")
    tx2gene = tmp_path / "tx2gene.tsv"
    tx2gene.write_text("transcript_id\tgene_id\tgene_name\nTX1\tG1\tG1\nTX2\tG1\tG1\n")
    output = tmp_path / "weighted.fastq"
    report = tmp_path / "weighting.txt"
    run_script("tile_transcripts.py", "--transcriptome", fasta, "--read-length", "500",
               "--step", "250", "--isoform-usage", usage, "--tx2gene", tx2gene,
               "--max-replicates", "3", "--weighting-report", report, "--output", output)
    names = [line[1:].strip() for line in output.read_text().splitlines() if line.startswith("@")]
    dominant = sum(1 for name in names if name.startswith("TX1_hap1"))
    minor = sum(1 for name in names if name.startswith("TX2_hap1"))
    assert dominant == 3 * minor
    assert report.read_text().startswith("expression")


def test_ase_input_keeps_libraries_without_reads(tmp_path: Path) -> None:
    """A library with no read for a gene measured zero and must stay in the fit."""
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample\tcondition\tdna_id\trna\nS1\tA\tD1\tx.bam\nS2\tB\tD2\ty.bam\n")
    counts = tmp_path / "counts.tsv"
    counts.write_text("sample\tID\tH1\tH2\tNonHS\nS1\tG1\t5\t3\t2\nS1\tG2\t1\t1\t1\nS2\tG2\t4\t4\t4\n")
    copy_number = tmp_path / "cn.tsv"
    copy_number.write_text("sample\tID\tCN_H1\tCN_H2\nD1\tG1\t1\t1\nD2\tG1\t2\t1\nD1\tG2\t1\t1\nD2\tG2\t1\t1\n")
    priors = tmp_path / "priors.tsv"
    priors.write_text("ID\tH1_prior\tH2_prior\nG1\t0.8\t0.8\nG2\t0.5\t0.5\n")
    output = tmp_path / "ase_input.tsv"
    run_script("prepare_ase_input.py", "--rna-counts", counts, "--samplesheet", sheet,
               "--copy-number", copy_number, "--priors", priors, "--output", output)
    rows = {(row["sample"], row["ID"]): row for row in read_tsv(output)}
    assert len(rows) == 4
    filled = rows[("S2", "G1")]
    assert (filled["H1_counts"], filled["H2_counts"], filled["NonHS_counts"]) == ("0", "0", "0")
    assert (filled["CN_H1"], filled["CN_H2"], filled["group"]) == ("2.0", "1.0", "B")

    # A combined table has copy number only on its rows: the missing library
    # takes it from another library of the same DNA sample.
    sheet.write_text("sample\tcondition\tdna_id\trna\nS1\tA\tD1\tx.bam\nS2\tA\tD1\ty.bam\n")
    combined = tmp_path / "combined.tsv"
    combined.write_text("sample\tID\tH1\tH2\tNonHS\tCN_H1\tCN_H2\nS1\tG1\t5\t3\t2\t3\t1\nS2\tG2\t4\t4\t4\t1\t1\n")
    run_script("prepare_ase_input.py", "--counts", combined, "--samplesheet", sheet,
               "--priors", priors, "--output", output)
    rows = {(row["sample"], row["ID"]): row for row in read_tsv(output)}
    assert (rows[("S2", "G1")]["CN_H1"], rows[("S2", "G1")]["H1_counts"]) == ("3.0", "0")


def test_split_ase_input_multiple_contrasts(tmp_path: Path) -> None:
    header = ("sample\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts"
              "\tCN_H1\tCN_H2\tH1_prior\tH2_prior\n")
    rows = []
    for group in ("A", "B", "C"):
        for replicate in range(1, 4):
            for gene in ("G1", "G2"):
                rows.append(f"{group}{replicate}\t{group}\t{gene}\t20\t20\t10\t1\t1\t0.8\t0.8\n")
    ase_input = tmp_path / "ASE_input.tsv"
    ase_input.write_text(header + "".join(rows))
    shards = tmp_path / "shards"
    run_script("split_ase_input.py", "--input", ase_input, "--contrast", "A:C,B:C",
               "--shards", "2", "--output-dir", shards)
    manifest = read_tsv(shards / "manifest.tsv")
    assert {row["contrast"] for row in manifest} == {"A_VS_C", "B_VS_C"}
    files = sorted(path.name for path in shards.glob("*.shard_*.tsv"))
    assert files == ["A_VS_C.shard_0001.tsv", "A_VS_C.shard_0002.tsv",
                     "B_VS_C.shard_0001.tsv", "B_VS_C.shard_0002.tsv"]
    # ASE.R refuses a shard holding more than one test, so each must be pure.
    for path in shards.glob("*.shard_*.tsv"):
        tests = {row["test"] for row in read_tsv(path)}
        groups = {row["group"] for row in read_tsv(path)}
        assert len(tests) == 1
        assert groups == set(tests.pop().split("_VS_"))


def test_split_ase_input_rejects_bad_contrasts(tmp_path: Path) -> None:
    header = ("sample\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts"
              "\tCN_H1\tCN_H2\tH1_prior\tH2_prior\n")
    body = "".join(f"A{i}\tA\tG1\t20\t20\t10\t1\t1\t0.8\t0.8\nB{i}\tB\tG1\t20\t20\t10\t1\t1\t0.8\t0.8\n"
                   for i in range(1, 4))
    ase_input = tmp_path / "ASE_input.tsv"
    ase_input.write_text(header + body)
    for contrast in ("A:B,A:B", "A:A", "A:B,A:Z", "A:B:C"):
        with pytest.raises(subprocess.CalledProcessError):
            run_script("split_ase_input.py", "--input", ase_input, "--contrast", contrast,
                       "--shards", "2", "--output-dir", tmp_path / f"out_{abs(hash(contrast))}")


def test_split_ase_input_rejects_pooled_sides(tmp_path: Path) -> None:
    """A contrast side is one condition: 'WT1+WT2' is not a condition."""
    header = ("sample\tgroup\tID\tH1_counts\tH2_counts\tNonHS_counts"
              "\tCN_H1\tCN_H2\tH1_prior\tH2_prior\n")
    rows = []
    for clone in ("WT1", "WT2", "GAIN"):
        for replicate in (1, 2, 3):
            rows.append(f"{clone}_{replicate}\t{clone}\tG1\t20\t20\t10\t1\t1\t0.8\t0.8\n")
    ase_input = tmp_path / "ASE_input.tsv"
    ase_input.write_text(header + "".join(rows))
    with pytest.raises(subprocess.CalledProcessError):
        run_script("split_ase_input.py", "--input", ase_input,
                   "--contrast", "GAIN:WT1+WT2", "--shards", "1",
                   "--output-dir", tmp_path / "bad")


def _copy_number_inputs(tmp_path: Path, chromosomes: dict, ploidy: str = "") -> dict:
    annotation = ["GID\tChrom\tSet"]
    hs, nonhs = [], []
    for chromosome, (h1, h2, shared) in chromosomes.items():
        for index in range(60):
            gene = f"{chromosome}_G{index:03d}"
            annotation.append(f"{gene}\t{chromosome}\tSingleton")
            hs += [f"{gene}_hap1\t{h1}", f"{gene}_hap2\t{h2}"]
            nonhs += [f"{gene}_hap1\t{shared}", f"{gene}_hap2\t{shared}"]
    (tmp_path / "ann.tsv").write_text("\n".join(annotation) + "\n")
    (tmp_path / "D1.HS.coverage.tsv").write_text("\n".join(hs) + "\n")
    (tmp_path / "D1.NonHS.coverage.tsv").write_text("\n".join(nonhs) + "\n")
    for group in ("HS_copy", "NonHS_copy"):
        (tmp_path / f"D1.{group}.coverage.tsv").write_text("")
    (tmp_path / "sheet.tsv").write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\n" f"S1\tA\tD1\t\t\t{ploidy}\n")
    return {path.name: path for path in tmp_path.glob("D1.*.coverage.tsv")}


def test_copy_number_calls_gains_and_losses(tmp_path: Path) -> None:
    chromosomes = {
        "chr1": (20, 20, 4), "chr2": (20, 20, 4), "chr3": (20, 20, 4),  # diploid
        "chr7": (40, 20, 4),                                            # gain, 2:1
        "chr9": (20, 0, 2),                                             # loss of H2
    }
    coverage = _copy_number_inputs(tmp_path, chromosomes, ploidy="chr7:2:1;chr9:1:0")
    run_script("compute_copy_number.py", "--samplesheet", tmp_path / "sheet.tsv",
               "--annotation-summary", tmp_path / "ann.tsv",
               "--coverage", *sorted(coverage.values()),
               "--output", tmp_path / "cn.tsv", "--qc-output", tmp_path / "qc.tsv",
               "--cn-change-threshold", "1.2")
    qc = {row["Chrom"]: row for row in read_tsv(tmp_path / "qc.tsv")}
    # Euploid chromosomes keep the neutral value whatever their allelic split.
    for chromosome in ("chr1", "chr2", "chr3"):
        assert qc[chromosome]["direction"] == "none"
        assert float(qc[chromosome]["CN_H1"]) == 1.0 == float(qc[chromosome]["CN_H2"])
    assert qc["chr7"]["direction"] == "gain"
    assert abs(float(qc["chr7"]["CN_H1"]) / float(qc["chr7"]["CN_H2"]) - 2) < 0.05
    # A loss is a departure too, and the missing haplotype never reaches zero.
    assert qc["chr9"]["direction"] == "loss"
    assert float(qc["chr9"]["CN_H2"]) > 0
    assert float(qc["chr9"]["CN_H2"]) < 0.1
    # The declaration only keeps a chromosome out of the baseline: chr9 was
    # declared 1:0 and still comes out of the depth, not out of the declaration.
    assert int(qc["chr1"]["baseline_genes"]) == 180
    assert abs(float(qc["chr9"]["Chrom_CN"]) - 1.0) < 0.05


def test_copy_number_gene_level_fallbacks(tmp_path: Path) -> None:
    """A gene without enough depth of its own leans on its chromosome."""
    annotation = ["GID\tChrom\tSet"]
    hs, nonhs = [], []
    for index in range(100):                  # chr1 is diploid and sets the baseline
        gene = f"chr1_G{index:03d}"
        annotation.append(f"{gene}\tchr1\tSingleton")
        hs += [f"{gene}_hap1\t20", f"{gene}_hap2\t20"]
        nonhs += [f"{gene}_hap1\t0", f"{gene}_hap2\t0"]
    for index in range(40):                   # chr2 carries three copies, split 0.75
        gene = f"chr2_G{index:03d}"
        annotation.append(f"{gene}\tchr2\tSingleton")
        hs += [f"{gene}_hap1\t45", f"{gene}_hap2\t15"]
        nonhs += [f"{gene}_hap1\t0", f"{gene}_hap2\t0"]
    # Two genes on that gained chromosome: one with no coverage at all, and one
    # with ample total depth but too few haplotype-specific reads to split it.
    annotation += ["NO_READS\tchr2\tComplex", "NO_SPLIT\tchr2\tComplex"]
    hs += ["NO_READS_hap1\t0", "NO_READS_hap2\t0", "NO_SPLIT_hap1\t1", "NO_SPLIT_hap2\t1"]
    nonhs += ["NO_READS_hap1\t0", "NO_READS_hap2\t0", "NO_SPLIT_hap1\t58", "NO_SPLIT_hap2\t58"]
    (tmp_path / "ann.tsv").write_text("\n".join(annotation) + "\n")
    (tmp_path / "D1.HS.coverage.tsv").write_text("\n".join(hs) + "\n")
    (tmp_path / "D1.NonHS.coverage.tsv").write_text("\n".join(nonhs) + "\n")
    for group in ("HS_copy", "NonHS_copy"):
        (tmp_path / f"D1.{group}.coverage.tsv").write_text("")
    (tmp_path / "sheet.tsv").write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\nS1\tA\tD1\t\t\t\n")

    run_script("compute_copy_number.py", "--samplesheet", tmp_path / "sheet.tsv",
               "--annotation-summary", tmp_path / "ann.tsv",
               "--coverage", *sorted(tmp_path.glob("D1.*.coverage.tsv")),
               "--output", tmp_path / "cn.tsv", "--qc-output", tmp_path / "qc.tsv",
               "--cn-resolution", "gene",
               "--min-haplotype-depth", "10", "--min-gene-depth", "1")
    cn = {row["ID"]: row for row in read_tsv(tmp_path / "cn.tsv")}
    qc = {row["Chrom"]: row for row in read_tsv(tmp_path / "qc.tsv")}
    assert qc["chr2"]["direction"] == "gain"
    assert abs(float(qc["chr2"]["Chrom_CN"]) - 3.0) < 0.1
    assert abs(float(qc["chr2"]["H1_prop"]) - 0.75) < 0.02

    # No coverage of any kind: the chromosome's per-haplotype ploidy, rather
    # than a homozygous deletion.
    assert float(cn["NO_READS"]["CN_H1"]) == float(qc["chr2"]["CN_H1"]) > 1.0
    assert float(cn["NO_READS"]["CN_H2"]) == float(qc["chr2"]["CN_H2"]) > 0.0

    # Ample total depth but no haplotype-specific reads: its own total ploidy,
    # split with the proportion seen in the singleton genes of its chromosome.
    total = float(cn["NO_SPLIT"]["CN_H1"]) + float(cn["NO_SPLIT"]["CN_H2"])
    assert abs(total - 3.0) < 0.1
    assert abs(float(cn["NO_SPLIT"]["CN_H1"]) / total - 0.75) < 0.02

    # A gene on the euploid chromosome is still measured on its own.
    assert float(cn["chr1_G000"]["CN_H1"]) == 1.0


def test_copy_number_gene_copies_option(tmp_path: Path) -> None:
    annotation = ["GID\tChrom\tSet"]
    hs, hs_copy = [], []
    for index in range(60):
        gene = f"chr1_G{index:03d}"
        annotation.append(f"{gene}\tchr1\tSingleton")
        hs += [f"{gene}_hap1\t20", f"{gene}_hap2\t20"]
        hs_copy += [f"{gene}_hap1\t0", f"{gene}_hap2\t0"]
    annotation.append("DUP\tchr1\tComplex")
    hs += ["DUP_hap1\t20", "DUP_hap2\t20"]
    hs_copy += ["DUP_1_hap1\t20", "DUP_1_hap2\t20"]        # a second copy elsewhere
    (tmp_path / "ann.tsv").write_text("\n".join(annotation) + "\n")
    (tmp_path / "D1.HS.coverage.tsv").write_text("\n".join(hs) + "\n")
    (tmp_path / "D1.HS_copy.coverage.tsv").write_text("\n".join(hs_copy) + "\n")
    for group in ("NonHS", "NonHS_copy"):
        (tmp_path / f"D1.{group}.coverage.tsv").write_text("")
    (tmp_path / "sheet.tsv").write_text(
        "sample\tcondition\tdna_id\trna\tdna\tploidy\nS1\tA\tD1\t\t\t\n")

    values = {}
    for mode in ("all", "primary"):
        run_script("compute_copy_number.py", "--samplesheet", tmp_path / "sheet.tsv",
                   "--annotation-summary", tmp_path / "ann.tsv",
                   "--coverage", *sorted(tmp_path.glob("D1.*.coverage.tsv")),
                   "--output", tmp_path / f"cn_{mode}.tsv",
                   "--qc-output", tmp_path / f"qc_{mode}.tsv",
                   "--cn-resolution", "gene", "--gene-copies", mode,
                   "--min-haplotype-depth", "10", "--min-gene-depth", "1")
        row = next(r for r in read_tsv(tmp_path / f"cn_{mode}.tsv") if r["ID"] == "DUP")
        values[mode] = float(row["CN_H1"]) + float(row["CN_H2"])
    # Summing the second copy doubles the dosage the model is told about.
    assert abs(values["all"] / values["primary"] - 2) < 0.05
