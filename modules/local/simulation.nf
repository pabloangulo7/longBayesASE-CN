process TILE_TRANSCRIPTS {
    tag "${meta.sample}"
    label 'process_low'

    input:
    tuple val(meta), val(hap), path(transcriptome)
    path transcript_quant
    path tx2gene
    val read_length
    val step
    val max_replicates

    output:
    tuple val(meta), path("${meta.sample}.fastq.gz"), emit: reads
    path 'weighting.txt', emit: weighting

    script:
    """
    tile_transcripts.py \
        --transcriptome ${transcriptome} \
        --read-length ${read_length} \
        --step ${step} \
        --transcript-quant ${transcript_quant} \
        --tx2gene ${tx2gene} \
        --max-replicates ${max_replicates} \
        --weighting-report weighting.txt \
        --output ${meta.sample}.fastq.gz
    """

    stub:
    """
    printf '@sim_read\nACGTACGT\n+\nIIIIIIII\n' | gzip -c > ${meta.sample}.fastq.gz
    echo uniform > weighting.txt
    """
}

process BUILD_MAPPING_PRIORS {
    tag 'mapping-bias priors'
    label 'process_low'

    input:
    path hap1_genes
    path hap2_genes
    path hap1_transcripts
    path hap2_transcripts
    path weighting
    val min_simulated_reads

    output:
    path 'mapping_priors.gene.tsv', emit: gene_priors
    path 'mapping_priors.transcript.tsv', emit: transcript_priors

    script:
    """
    weighting=\$(cat ${weighting})
    build_priors.py \
        --hap1-counts ${hap1_genes} \
        --hap2-counts ${hap2_genes} \
        --min-simulated-reads ${min_simulated_reads} \
        --weighting "\${weighting}" \
        --output mapping_priors.gene.tsv
    build_priors.py \
        --hap1-counts ${hap1_transcripts} \
        --hap2-counts ${hap2_transcripts} \
        --min-simulated-reads ${min_simulated_reads} \
        --weighting "\${weighting}" \
        --output mapping_priors.transcript.tsv
    """

    stub:
    """
    printf '# weighting: uniform\nID\tH1_prior\tH2_prior\tH1_misassignment\tH2_misassignment\tn_hap1\tn_hap2\nGENE1\t0.8\t0.8\t0\t0\t100\t100\n' > mapping_priors.gene.tsv
    cp mapping_priors.gene.tsv mapping_priors.transcript.tsv
    """
}
