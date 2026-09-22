process DNA_ALIGN {
    tag "${meta.dna_id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path diploid_fasta

    output:
    tuple val(meta), path("${meta.dna_id}.dna.bam"), path("${meta.dna_id}.dna.bam.bai"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.dna_id}\\tPU:Unknown\\tPL:ONT\\tLB:${meta.dna_id}\\tSM:${meta.dna_id}"
    """
    minimap2 \
        -y -ax lr:hq \
        -N 200 \
        -R '${rg}' \
        -t ${task.cpus} \
        ${diploid_fasta} ${reads} \
      | samtools sort -@ ${task.cpus} -o ${meta.dna_id}.dna.bam -
    samtools index -@ ${task.cpus} ${meta.dna_id}.dna.bam
    """

    stub:
    """
    touch ${meta.dna_id}.dna.bam ${meta.dna_id}.dna.bam.bai
    """
}

process DNA_MAX_SCORE {
    tag "${meta.dna_id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.dna_id}.maxscore.bam"), path("${meta.dna_id}.maxscore.bam.bai"), path("${meta.dna_id}.alignment_readgroups.tsv"), emit: maxscore

    script:
    """
    # The classifier only needs each read's alignments to be adjacent, which
    # collate gives without a full name sort.
    samtools collate -@ ${task.cpus} -o input.namesort.bam ${bam}
    classify_dna_alignments.py \
        --bam input.namesort.bam \
        --output-bam ${meta.dna_id}.maxscore.namesort.bam \
        --readgroups ${meta.dna_id}.alignment_readgroups.tsv \
        --threads ${task.cpus}
    samtools sort -@ ${task.cpus} -o ${meta.dna_id}.maxscore.bam ${meta.dna_id}.maxscore.namesort.bam
    samtools index -@ ${task.cpus} ${meta.dna_id}.maxscore.bam
    rm input.namesort.bam ${meta.dna_id}.maxscore.namesort.bam
    """

    stub:
    """
    touch ${meta.dna_id}.maxscore.bam ${meta.dna_id}.maxscore.bam.bai
    printf 'Read_ID\tH1\tH2\tGroup\tMax_Score\tUnclassified\n' > ${meta.dna_id}.alignment_readgroups.tsv
    """
}

process DNA_FEATURE_ASSIGN {
    tag "${meta.dna_id}"
    label 'process_high'

    input:
    tuple val(meta), path(coord_bam, stageAs: 'maxscore_input/*'), path(bai), path(alignment_groups)
    path genes_saf

    output:
    tuple val(meta), path("${meta.dna_id}.maxscore.bam"), path("${meta.dna_id}.readgroups_genes.tsv"), emit: split_input
    tuple val(meta), path("${meta.dna_id}.gene_counts.tsv"), emit: counts
    tuple val(meta), path("${meta.dna_id}.mapping_qc.tsv"), path("${meta.dna_id}.complex_reads.tsv"), emit: qc

    script:
    """
    featureCounts \
        -T ${task.cpus} -O -Q 0 -L -M -R CORE -F SAF -s 0 \
        -a ${genes_saf} \
        -o ${meta.dna_id}.featureCounts.txt \
        ${coord_bam}
    classify_dna_genes.py \
        --featurecounts ${meta.dna_id}.maxscore.bam.featureCounts \
        --readgroups ${alignment_groups} \
        --output-prefix ${meta.dna_id}
    ln -L ${coord_bam} ${meta.dna_id}.maxscore.bam
    """

    stub:
    """
    touch ${meta.dna_id}.maxscore.bam
    printf 'Read_ID\tGroup\n' > ${meta.dna_id}.readgroups_genes.tsv
    printf 'GID\tH1\tH1_multigene\tH1_multimapping\tH1_multimapping_multigene\tH1_multimapping_multigene_copy\tH1_multimapping_multigene_diff\tH2\tH2_multigene\tH2_multimapping\tH2_multimapping_multigene\tH2_multimapping_multigene_copy\tH2_multimapping_multigene_diff\tNonHS\tNonHS_multigene\tNonHS_multimapping\tNonHS_multimapping_multigene\tNonHS_multigene_copy\tNonHS_multimapping_multigene_copy\tNonHS_multimapping_multigene_diff\tNonHS_complex\tNonHS_multimapping_complex\tUnclassified\nGENE1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' > ${meta.dna_id}.gene_counts.tsv
    printf 'category\treads\n' > ${meta.dna_id}.mapping_qc.tsv
    printf 'Read_ID\tGroup\tGenes\n' > ${meta.dna_id}.complex_reads.tsv
    """
}

process DNA_SPLIT {
    tag "${meta.dna_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(coord_bam), path(readgroups)

    output:
    tuple val(meta), path("${meta.dna_id}.HS.bam"), path("${meta.dna_id}.NonHS.bam"), path("${meta.dna_id}.HS_copy.bam"), path("${meta.dna_id}.NonHS_copy.bam"), emit: bams

    script:
    """
    split_dna_bam.py \
        --bam ${coord_bam} \
        --readgroups ${readgroups} \
        --output-prefix split \
        --threads ${task.cpus}
    for group in HS NonHS HS_copy NonHS_copy; do
        samtools sort -@ ${task.cpus} -o ${meta.dna_id}.\${group}.bam split.\${group}.bam
    done
    """

    stub:
    """
    touch ${meta.dna_id}.HS.bam ${meta.dna_id}.NonHS.bam ${meta.dna_id}.HS_copy.bam ${meta.dna_id}.NonHS_copy.bam
    """
}

process DNA_COVERAGE {
    tag "${meta.dna_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(hs_bam), path(nonhs_bam), path(hs_copy_bam), path(nonhs_copy_bam)
    path genes_bed

    output:
    tuple val(meta), path("${meta.dna_id}.*.coverage.tsv"), emit: coverage

    script:
    """
    # The four read groups are independent, so they run in parallel.
    for group in HS NonHS HS_copy NonHS_copy; do
        bedtools coverage -a ${genes_bed} -b ${meta.dna_id}.\${group}.bam -mean \
          | awk 'BEGIN{OFS="\\t"}{print \$4,\$NF}' > ${meta.dna_id}.\${group}.coverage.tsv &
    done
    wait
    """

    stub:
    """
    for group in HS NonHS HS_copy NonHS_copy; do
        printf 'GENE1_hap1\t10\nGENE1_hap2\t10\n' > ${meta.dna_id}.\${group}.coverage.tsv
    done
    """
}

process COMPUTE_COPY_NUMBER {
    tag 'copy number'
    label 'process_medium'

    input:
    path samplesheet
    path annotation_summary
    path coverage_files
    val min_haplotype_depth
    val min_gene_depth
    val gene_copies
    val cn_change_threshold
    val cn_resolution

    output:
    path 'gene_copy_numbers.tsv', emit: copy_number
    path 'copy_number.qc.tsv', emit: qc

    script:
    """
    compute_copy_number.py \
        --samplesheet ${samplesheet} \
        --annotation-summary ${annotation_summary} \
        --coverage ${coverage_files} \
        --min-haplotype-depth ${min_haplotype_depth} \
        --min-gene-depth ${min_gene_depth} \
        --gene-copies ${gene_copies} \
        --cn-change-threshold ${cn_change_threshold} \
        --cn-resolution '${cn_resolution}' \
        --output gene_copy_numbers.tsv \
        --qc-output copy_number.qc.tsv
    """

    stub:
    """
    printf 'sample\tID\tChrom\tCN_H1\tCN_H2\nDNA1\tGENE1\tchr1\t1\t1\n' > gene_copy_numbers.tsv
    printf 'sample\tChrom\tbaseline_depth\tbaseline_genes\tchromosome_singleton_genes\thaplotype_informative_genes\texpected_ploidy\tChrom_CN\tratio\tchanged\tdirection\tH1_prop\tH2_prop\tCN_H1\tCN_H2\nDNA1\tchr1\t10\t1\t1\t1\tFalse\t2\t1.0\tFalse\tnone\t0.5\t0.5\t1\t1\n' > copy_number.qc.tsv
    """
}

process COPY_NUMBER_FROM_PLOIDY {
    tag 'copy number from expected ploidy'
    label 'process_low'

    input:
    path samplesheet
    path annotation_summary

    output:
    path 'gene_copy_numbers.tsv', emit: copy_number

    script:
    """
    copy_number_from_ploidy.py \
        --samplesheet ${samplesheet} \
        --annotation-summary ${annotation_summary} \
        --output gene_copy_numbers.tsv
    """

    stub:
    """
    printf 'sample\tID\tChrom\tCN_H1\tCN_H2\nDNA1\tGENE1\tchr1\t1\t1\n' > gene_copy_numbers.tsv
    """
}
