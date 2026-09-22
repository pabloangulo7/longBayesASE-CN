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

process OARFISH_ASSIGN {
    tag "${meta.sample}"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    path tx2gene
    val score_threshold
    val display_threshold
    val strand

    output:
    tuple val(meta), path("${meta.sample}.assignments.tsv"), emit: assignments
    tuple val(meta), path("${meta.sample}.oarfish.tsv"), emit: quantification

    script:
    def strandFlag = strand ? "--strand-filter '${strand}'" : ''
    """
    oarfish \
        --verbose \
        -j ${task.cpus} \
        --alignments ${bam} \
        --output ${meta.sample}.oarfish \
        --filter-group no-filters \
        ${strandFlag} \
        --score-threshold ${score_threshold} \
        --display-thresh ${display_threshold} \
        --model-coverage \
        --write-assignment-probs=uncompressed
    if [ -s ${meta.sample}.oarfish.prob ]; then
        mv ${meta.sample}.oarfish.prob ${meta.sample}.prob
    elif [ ! -s ${meta.sample}.prob ]; then
        echo 'Oarfish did not produce an assignment-probability file' >&2
        exit 1
    fi
    if [ ! -s ${meta.sample}.oarfish.quant ]; then
        echo 'Oarfish did not produce its quantification table' >&2
        exit 1
    fi
    mv ${meta.sample}.oarfish.quant ${meta.sample}.oarfish.tsv
    parse_oarfish_prob.py \
        --prob ${meta.sample}.prob \
        --tx2gene ${tx2gene} \
        --output ${meta.sample}.assignments.tsv
    """

    stub:
    """
    printf 'read_id\ttranscript_id\tgene_id\tprobability\nstub_read\tTX1_hap1\tGENE1_hap1\t1\n' > ${meta.sample}.assignments.tsv
    printf 'transcript\test_counts\nTX1_hap1\t1\n' > ${meta.sample}.oarfish.tsv
    """
}

process RNA_HAPLOTYPE_COUNT {
    tag "${meta.sample}"
    label 'process_medium'

    input:
    tuple val(meta), path(assignments)
    val gene_resolve_probability
    val isoform_resolve_probability
    val gene_min_probability
    val isoform_min_probability

    output:
    tuple val(meta), path("${meta.sample}.genes.counts.tsv"), emit: gene_counts
    tuple val(meta), path("${meta.sample}.isoforms.counts.tsv"), emit: isoform_counts
    tuple val(meta), path("${meta.sample}.*.qc.tsv"), path("${meta.sample}.*.complex_reads.tsv"), emit: qc
    tuple val(meta), path("${meta.sample}.genes.readgroups.tsv"), path("${meta.sample}.isoforms.readgroups.tsv"), emit: readgroups

    script:
    """
    # Gene and isoform level are independent passes over the same file. Each one
    # is waited for by PID, so a pass killed for memory fails the task and is retried.
    pids=()
    for mode in genes isoforms; do
        if [ "\${mode}" = genes ]; then
            resolve=${gene_resolve_probability}
            min_probability=${gene_min_probability}
        else
            resolve=${isoform_resolve_probability}
            min_probability=${isoform_min_probability}
        fi
        count_rna_haplotypes.py \
            --assignments ${assignments} \
            --mode \${mode} \
            --resolve-threshold \${resolve} \
            --min-probability \${min_probability} \
            --output-prefix ${meta.sample}.\${mode} &
        pids+=(\$!)
    done
    for pid in "\${pids[@]}"; do wait "\${pid}"; done
    """

    stub:
    """
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nGENE1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.genes.counts.tsv
    printf 'ID\tH1\tH1_multimapping\tH2\tH2_multimapping\tNonHS\tNonHS_multimapping\nTX1\t10\t0\t10\t0\t5\t0\n' > ${meta.sample}.isoforms.counts.tsv
    printf 'metric\treads\nresolved\t25\n' > ${meta.sample}.genes.qc.tsv
    cp ${meta.sample}.genes.qc.tsv ${meta.sample}.isoforms.qc.tsv
    printf 'Read_ID\tGroup\tResolved\n' > ${meta.sample}.genes.readgroups.tsv
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
