process RNA_ALIGN {
    tag "${meta.sample}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path transcriptome

    output:
    tuple val(meta), path("${meta.sample}.rna.bam"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.sample}\\tPU:Unknown\\tPL:ONT\\tLB:${meta.sample}\\tSM:${meta.sample}"
    """
    minimap2 \
        -y --eqx -ax lr:hq \
        -N 200 \
        -R '${rg}' \
        -t ${task.cpus} \
        ${transcriptome} ${reads} \
      | samtools view -@ ${task.cpus} -b - \
      > ${meta.sample}.rna.bam
    """

    stub:
    """
    touch ${meta.sample}.rna.bam
    """
}

process OARFISH_QUANT {
    tag "${meta.sample}"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    val score_threshold
    val strand
    val bootstraps

    output:
    tuple val(meta), path("${meta.sample}.quant"), emit: quant
    tuple val(meta), path("${meta.sample}.ambig_info.tsv"), path("${meta.sample}.meta_info.json"), emit: info
    tuple val(meta), path("${meta.sample}.infreps.pq"), emit: infreps, optional: true

    script:
    """
    oarfish \
        -j ${task.cpus} \
        --alignments ${bam} \
        --output ${meta.sample} \
        --filter-group no-filters \
        --strand-filter ${strand} \
        --score-threshold ${score_threshold} \
        --model-coverage \
        --num-bootstraps ${bootstraps}
    """

    stub:
    """
    printf 'tname\tlen\tnum_reads\nTX1_hap1\t1000\t15\nTX1_hap2\t1000\t10\n' > ${meta.sample}.quant
    printf 'unique_reads\tambig_reads\ttotal_reads\n10\t5\t15\n5\t5\t10\n' > ${meta.sample}.ambig_info.tsv
    echo '{}' > ${meta.sample}.meta_info.json
    """
}

process RNA_HAPLOTYPE_COUNT {
    tag "${meta.sample}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)
    path tx2gene
    val strand

    output:
    tuple val(meta), path("${meta.sample}.genes.counts.tsv"), emit: gene_counts
    tuple val(meta), path("${meta.sample}.isoforms.counts.tsv"), emit: isoform_counts
    tuple val(meta), path("${meta.sample}.*.qc.tsv"), path("${meta.sample}.*.complex_reads.tsv"), emit: qc
    tuple val(meta), path("${meta.sample}.genes.readgroups.tsv"), path("${meta.sample}.isoforms.readgroups.tsv"), emit: readgroups

    script:
    """
    count_rna_haplotypes.py \
        --bam ${bam} \
        --tx2gene ${tx2gene} \
        --strand ${strand} \
        --threads ${task.cpus} \
        --output-prefix ${meta.sample}
    """

    stub:
    """
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nGENE1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.genes.counts.tsv
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nTX1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.isoforms.counts.tsv
    printf 'metric\treads\nH1\t10\n' > ${meta.sample}.genes.qc.tsv
    cp ${meta.sample}.genes.qc.tsv ${meta.sample}.isoforms.qc.tsv
    printf 'Read_ID\tGroup\n' > ${meta.sample}.genes.readgroups.tsv
    cp ${meta.sample}.genes.readgroups.tsv ${meta.sample}.isoforms.readgroups.tsv
    touch ${meta.sample}.genes.complex_reads.tsv ${meta.sample}.isoforms.complex_reads.tsv
    """
}

process MERGE_RNA_COUNTS {
    tag 'unified RNA counts'
    label 'process_low'

    input:
    path counts
    val outname

    output:
    path outname, emit: counts

    script:
    """
    merge_rna_counts.py --counts ${counts} --output ${outname}
    """

    stub:
    """
    printf 'sample\tID\tH1\tH2\tNonHS\nSample1\tGENE1\t10\t10\t5\n' > ${outname}
    """
}

process MERGE_OARFISH_QUANT {
    tag 'Oarfish expression tables'
    label 'process_low'

    input:
    path quant
    path tx2gene

    output:
    path 'gene_quant.tsv', emit: genes
    path 'transcript_quant.tsv', emit: transcripts

    script:
    """
    merge_oarfish_quant.py \
        --quant ${quant} \
        --tx2gene ${tx2gene} \
        --genes-output gene_quant.tsv \
        --transcripts-output transcript_quant.tsv
    """

    stub:
    """
    printf 'sample\tID\tnum_reads\nSample1\tGENE1\t25\n' > gene_quant.tsv
    printf 'sample\tID\tnum_reads\nSample1\tTX1\t25\n' > transcript_quant.tsv
    """
}
