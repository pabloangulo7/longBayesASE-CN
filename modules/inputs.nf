process VALIDATE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_low'

    input:
    path samplesheet
    val step

    output:
    path 'samplesheet.validated.tsv', emit: samplesheet

    script:
    """
    validate_samplesheet.py \
        --input ${samplesheet} \
        --output samplesheet.validated.tsv \
        --step ${step}
    """

    stub:
    """
    cp ${samplesheet} samplesheet.validated.tsv
    """
}

process NORMALIZE_READS {
    tag "${meta.sample ?: meta.dna_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path('reads.fastq.gz'), emit: reads

    script:
    """
    : > reads.fastq.gz
    for input_read in input/*; do
      case "\${input_read}" in
        *.bam|*.ubam|*.cram)
          samtools fastq -@ ${task.cpus} -T '*' "\${input_read}" | gzip -c >> reads.fastq.gz
          ;;
        *.fastq.gz|*.fq.gz)
          cat "\${input_read}" >> reads.fastq.gz
          ;;
        *.fastq|*.fq)
          gzip -c "\${input_read}" >> reads.fastq.gz
          ;;
        *)
          echo "Unsupported read format: \${input_read}" >&2
          exit 1
          ;;
      esac
    done
    """

    stub:
    """
    printf '@stub_read\nACGTACGT\n+\nFFFFFFFF\n' | gzip -c > reads.fastq.gz
    """
}

process PREPARE_RNA_BAM {
    tag "${meta.sample}"
    label 'process_low'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.rna.bam"), emit: bam

    script:
    """
    ln -L ${bam} ${meta.sample}.rna.bam
    """

    stub:
    """
    touch ${meta.sample}.rna.bam
    """
}
