process PREPARE_REFERENCE {
    tag 'diploid reference'
    label 'process_medium'

    input:
    path diploid_fasta

    output:
    path 'reference/hap1.fa', emit: hap1_fasta
    path 'reference/hap2.fa', emit: hap2_fasta
    path 'reference/hap1.chroms.csv', emit: hap1_map
    path 'reference/hap2.chroms.csv', emit: hap2_map
    path 'reference/reference_summary.tsv', emit: summary
    path 'reference/excluded_contigs.tsv', emit: excluded

    script:
    """
    prepare_reference.py --fasta ${diploid_fasta} --output-dir reference
    """

    stub:
    """
    mkdir -p reference
    printf '>chr1_hap1\nACGTACGTACGT\n' > reference/hap1.fa
    printf '>chr1_hap2\nACGTACGTACGT\n' > reference/hap2.fa
    printf 'chr1,chr1_hap1\n' > reference/hap1.chroms.csv
    printf 'chr1,chr1_hap2\n' > reference/hap2.chroms.csv
    printf 'haplotype\tcontigs\nhap1\t1\nhap2\t1\n' > reference/reference_summary.tsv
    printf 'contig\treason\n' > reference/excluded_contigs.tsv
    """
}

process SPLIT_DIPLOID_GFF {
    tag 'existing diploid annotation'
    label 'process_low'

    input:
    path gff3

    output:
    path 'hap1.gff3', emit: hap1
    path 'hap2.gff3', emit: hap2

    script:
    """
    split_diploid_gff.py --gff3 ${gff3} --hap1 hap1.gff3 --hap2 hap2.gff3
    """

    stub:
    """
    printf '##gff-version 3\nchr1_hap1\tstub\tgene\t1\t12\t.\t+\t.\tID=GENE1;gene_id=GENE1\nchr1_hap1\tstub\ttranscript\t1\t12\t.\t+\t.\tID=TX1;Parent=GENE1;gene_id=GENE1;transcript_id=TX1\nchr1_hap1\tstub\texon\t1\t12\t.\t+\t.\tID=EX1;Parent=TX1;gene_id=GENE1;transcript_id=TX1\n' > hap1.gff3
    sed 's/_hap1/_hap2/g' hap1.gff3 > hap2.gff3
    """
}

process LIFTOFF_ANNOTATE {
    tag "${hap}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(hap), path(target_fasta), path(chromosome_map)
    path source_fasta
    path source_annotation

    output:
    tuple val(hap), path("${hap}.gff3"), path("${hap}.unmapped.txt"), emit: annotation

    script:
    """
    liftoff \
        -p ${task.cpus} \
        -polish \
        -exclude_partial \
        -copies \
        -sc 0.90 \
        -chroms ${chromosome_map} \
        -g ${source_annotation} \
        -o ${hap}.gff3 \
        -u ${hap}.unmapped.txt \
        ${target_fasta} ${source_fasta}
    """

    stub:
    """
    printf '##gff-version 3\nchr1_${hap}\tstub\tgene\t1\t12\t.\t+\t.\tID=GENE1;gene_id=GENE1;gene_name=GENE1;extra_copy_number=0\nchr1_${hap}\tstub\ttranscript\t1\t12\t.\t+\t.\tID=TX1;Parent=GENE1;gene_id=GENE1;transcript_id=TX1;extra_copy_number=0\nchr1_${hap}\tstub\texon\t1\t12\t.\t+\t.\tID=EX1;Parent=TX1;gene_id=GENE1;transcript_id=TX1;extra_copy_number=0\n' > ${hap}.gff3
    touch ${hap}.unmapped.txt
    """
}

process BUILD_ANNOTATION_BUNDLE {
    tag 'annotation bundle'
    label 'process_medium'

    input:
    path hap1_gff
    path hap2_gff
    path diploid_fasta

    output:
    path 'bundle/annotation_summary.tsv', emit: summary
    path 'bundle/genes.bed', emit: bed
    path 'bundle/genes.saf', emit: saf
    path 'bundle/tx2gene.tsv', emit: tx2gene
    path 'bundle/diploid.gff3', emit: gff
    path 'bundle/diploid_transcriptome.fa', emit: transcriptome
    path 'bundle/hap1_transcriptome.fa', emit: hap1_transcriptome
    path 'bundle/hap2_transcriptome.fa', emit: hap2_transcriptome
    path 'bundle/transcriptome_summary.tsv', emit: transcriptome_summary

    script:
    """
    mkdir -p bundle
    annotation_bundle.py \
        --hap1-gff ${hap1_gff} \
        --hap2-gff ${hap2_gff} \
        --output-dir bundle
    gffread ${hap1_gff} -g ${diploid_fasta} -w hap1.raw.fa
    gffread ${hap2_gff} -g ${diploid_fasta} -w hap2.raw.fa
    filter_transcriptome.py --hap1 hap1.raw.fa --hap2 hap2.raw.fa --output-dir bundle
    """

    stub:
    """
    mkdir -p bundle
    printf 'GID\tgene_name\tChrom\tcopies_h1\tcopies_h2\tSet\nGENE1\tGENE1\tchr1\t1\t1\tSingleton\n' > bundle/annotation_summary.tsv
    printf 'chr1_hap1\t0\t12\tGENE1_hap1\nchr1_hap2\t0\t12\tGENE1_hap2\n' > bundle/genes.bed
    printf 'GeneID\tChr\tStart\tEnd\tStrand\nGENE1_hap1\tchr1_hap1\t1\t12\t+\nGENE1_hap2\tchr1_hap2\t1\t12\t+\n' > bundle/genes.saf
    printf 'transcript_id\tgene_id\tgene_name\nTX1\tGENE1\tGENE1\n' > bundle/tx2gene.tsv
    printf '##gff-version 3\n' > bundle/diploid.gff3
    printf '>TX1_hap1\nACGTACGTACGT\n>TX1_hap2\nACGTACGTACGT\n' > bundle/diploid_transcriptome.fa
    printf '>TX1_hap1\nACGTACGTACGT\n' > bundle/hap1_transcriptome.fa
    printf '>TX1_hap2\nACGTACGTACGT\n' > bundle/hap2_transcriptome.fa
    printf 'haplotype\ttranscripts\nhap1\t1\nhap2\t1\n' > bundle/transcriptome_summary.tsv
    """
}

process BUILD_FEATURE_MAP {
    tag 'feature map'
    label 'process_low'

    input:
    path gff3

    output:
    path 'map/gene_map.tsv', emit: genes
    path 'map/tx2gene.tsv', emit: tx2gene

    script:
    """
    build_feature_map.py --gff3 ${gff3} --output-dir map
    """

    stub:
    """
    mkdir -p map
    printf 'GID\tChrom\tgene_name\tSet\nGENE1\tchr1\tGENE1\tNA\n' > map/gene_map.tsv
    printf 'transcript_id\tgene_id\tgene_name\nTX1\tGENE1\tGENE1\n' > map/tx2gene.tsv
    """
}
