process SIMULATE_PRIOR_READS {
    tag "${meta.sample}"
    label 'process_low'

    input:
    tuple val(meta), path(transcriptome)
    path read_intervals

    output:
    tuple val(meta), path("${meta.sample}.fastq.gz"), emit: reads

    script:
    """
    simulate_prior_reads.py \
        --intervals ${read_intervals} \
        --transcriptome ${transcriptome} \
        --output ${meta.sample}.fastq.gz
    """

    stub:
    """
    printf '@sim_read\nACGTACGT\n+\nIIIIIIII\n' | gzip -c > ${meta.sample}.fastq.gz
    """
}

process BUILD_MAPPING_PRIORS {
    tag "${level}"
    label 'process_low'

    input:
    path hap1_counts
    path hap2_counts
    val level
    val min_reads

    output:
    path "mapping_priors.${level}.tsv", emit: priors

    script:
    """
    build_priors.py \
        --hap1-counts ${hap1_counts} \
        --hap2-counts ${hap2_counts} \
        --min-reads ${min_reads} \
        --output mapping_priors.${level}.tsv
    """

    stub:
    """
    printf 'ID\tH1_prior\tH2_prior\tH1_misassignment\tH2_misassignment\tn_hap1\tn_hap2\nFEATURE1\t0.8\t0.8\t0\t0\t100\t100\n' > mapping_priors.${level}.tsv
    """
}
