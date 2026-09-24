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
    val intervals_per_feature

    output:
    tuple val(meta), path("${meta.sample}.gene_HS_counts.tsv"), emit: gene_hs_counts
    tuple val(meta), path("${meta.sample}.transcript_HS_counts.tsv"), emit: transcript_hs_counts
    tuple val(meta), path("${meta.sample}.*_qc.tsv"), path("${meta.sample}.*_complex_reads.tsv"), emit: qc
    tuple val(meta), path("${meta.sample}.gene_readgroups.tsv"), path("${meta.sample}.transcript_readgroups.tsv"), emit: readgroups
    tuple val(meta), path("${meta.sample}.gene_read_intervals.tsv.gz"), emit: gene_read_intervals, optional: true
    tuple val(meta), path("${meta.sample}.transcript_read_intervals.tsv.gz"), emit: transcript_read_intervals, optional: true

    script:
    """
    count_rna_haplotypes.py \
        --bam ${bam} \
        --tx2gene ${tx2gene} \
        --strand ${strand} \
        --intervals-per-feature ${intervals_per_feature} \
        --threads ${task.cpus} \
        --output-prefix ${meta.sample}
    """

    stub:
    """
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nGENE1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.gene_HS_counts.tsv
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nTX1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.transcript_HS_counts.tsv
    printf 'metric\treads\nH1\t10\n' > ${meta.sample}.gene_qc.tsv
    cp ${meta.sample}.gene_qc.tsv ${meta.sample}.transcript_qc.tsv
    printf 'Read_ID\tGroup\n' > ${meta.sample}.gene_readgroups.tsv
    cp ${meta.sample}.gene_readgroups.tsv ${meta.sample}.transcript_readgroups.tsv
    touch ${meta.sample}.gene_complex_reads.tsv ${meta.sample}.transcript_complex_reads.tsv
    if [ ${intervals_per_feature} -gt 0 ]; then
        printf 'ID\ttranscript\tstart\tend\nGENE1\tTX1\t0\t8\n' | gzip -c > ${meta.sample}.gene_read_intervals.tsv.gz
        printf 'ID\ttranscript\tstart\tend\nTX1\tTX1\t0\t8\n' | gzip -c > ${meta.sample}.transcript_read_intervals.tsv.gz
    fi
    """
}

process MERGE_HS_COUNTS {
    tag "${outname}"
    label 'process_low'

    input:
    path counts
    val outname

    output:
    path outname, emit: counts

    script:
    """
    merge_hs_counts.py --counts ${counts} --output ${outname}
    """

    stub:
    """
    printf 'sample\tID\tH1\tH2\tNonHS\nSample1\tGENE1\t10\t10\t5\n' > ${outname}
    """
}

process MERGE_READ_INTERVALS {
    tag "${level}"
    label 'process_low'

    input:
    path intervals
    val level

    output:
    path "${level}_read_intervals.tsv.gz", emit: intervals

    script:
    """
    {
        printf 'sample\tID\ttranscript\tstart\tend\n'
        for file in \$(ls ${intervals} | sort); do
            sample=\${file%.${level}_read_intervals.tsv.gz}
            zcat "\${file}" | tail -n +2 | awk -v sample="\${sample}" 'BEGIN { OFS = "\t" } { print sample, \$0 }'
        done
    } | gzip -c > ${level}_read_intervals.tsv.gz
    """

    stub:
    """
    printf 'sample\tID\ttranscript\tstart\tend\nSample1\tFEATURE1\tTX1\t0\t8\n' | gzip -c > ${level}_read_intervals.tsv.gz
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
