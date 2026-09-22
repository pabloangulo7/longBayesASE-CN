#!/usr/bin/env nextflow

include { VALIDATE_SAMPLESHEET; NORMALIZE_READS as NORMALIZE_DNA_READS; NORMALIZE_READS as NORMALIZE_RNA_READS; PREPARE_DNA_BAM; PREPARE_RNA_BAM } from './modules/local/inputs'
include { PREPARE_REFERENCE; SPLIT_DIPLOID_GFF; LIFTOFF_ANNOTATE; BUILD_ANNOTATION_BUNDLE; BUILD_FEATURE_MAP } from './modules/local/annotation'
include { DNA_ALIGN; DNA_MAX_SCORE; DNA_FEATURE_ASSIGN; DNA_SPLIT; DNA_COVERAGE; COMPUTE_COPY_NUMBER; COPY_NUMBER_FROM_PLOIDY } from './modules/local/dna'
include { RNA_ALIGN; OARFISH_ASSIGN; RNA_HAPLOTYPE_COUNT; MERGE_RNA_COUNTS; MERGE_RNA_COUNTS as MERGE_ISOFORM_COUNTS } from './modules/local/rna'
include { RNA_ALIGN as SIM_RNA_ALIGN; OARFISH_ASSIGN as SIM_OARFISH_ASSIGN; RNA_HAPLOTYPE_COUNT as SIM_RNA_HAPLOTYPE_COUNT } from './modules/local/rna'
include { TILE_TRANSCRIPTS; BUILD_MAPPING_PRIORS } from './modules/local/simulation'
include { PREPARE_ASE_INPUT; PREPARE_COMBINED_ASE_INPUT; COMPILE_ASE_MODEL; SPLIT_ASE_INPUT; RUN_ASE_SHARD; MERGE_ASE_RESULTS } from './modules/local/ase'


def validateParams() {
    def schema = new groovy.json.JsonSlurper().parse(file("${projectDir}/nextflow_schema.json"))
    // Each line ends with its operator: Nextflow 26 rejects a line that opens with one.
    def known = ((schema.definitions ?: [:]).values().collectMany { definition -> (definition.properties ?: [:]).keySet() } +
                 (schema.properties ?: [:]).keySet() +
                 params.keySet().findAll { name -> name.contains('-') }) as Set
    def unknown = params.keySet().findAll { name -> !known.contains(name) && !name.contains('-') }
    if (unknown) error "Unknown parameter(s): ${unknown.sort().join(', ')}"
}

def requiredFile(value, label) {
    if (!value) error "Missing required parameter --${label}"
    channel.value(file(value, checkIfExists: true))
}

def optionalFile(value, fallback) {
    channel.value(file(value ?: fallback, checkIfExists: true))
}

def emptyTx2gene() { optionalFile(null, "${projectDir}/assets/empty_tx2gene.tsv") }

// Absent --isoform_usage resolves to a header-only table, so the tiler falls
// back to one tile per transcript and reproduces the unweighted priors exactly.
// At isoform level the weighting scales a transcript's numerator and denominator
// by the same factor and cancels, so it is never applied there.
def isoformUsage() {
    if (params.level == 'isoform' && params.isoform_usage) {
        log.warn 'Ignoring --isoform_usage: isoform-level priors do not average over isoforms, so the weighting cancels out'
    }
    optionalFile(params.level == 'isoform' ? null : params.isoform_usage,
                 "${projectDir}/assets/empty_isoform_usage.tsv")
}

def isoformLevel() { params.level == 'isoform' }

// Pick one haplotype's annotation out of the Liftoff channel.
def haplotypeAnnotation(channel, wanted) {
    channel
        .filter { hap, _gff, _unmapped -> hap == wanted }
        .map { _hap, gff, _unmapped -> gff }
}

// Pick one simulated library's count table out of the per-sample channel.
def simulatedCounts(channel, wanted) {
    channel
        .filter { meta, _path -> meta.sample == wanted }
        .map { _meta, path -> path }
}

def samplesHaveDNA(path) {
    def lines = file(path, checkIfExists: true).readLines()
    def header = lines[0].split('\t', -1) as List
    def index = header.indexOf('dna')
    return index >= 0 && lines.drop(1).any { line -> line.split('\t', -1)[index].trim() }
}


workflow REFERENCE_WF {
    take:
    fasta

    main:
    if (params.gff3) {
        SPLIT_DIPLOID_GFF(requiredFile(params.gff3, 'gff3'))
        BUILD_ANNOTATION_BUNDLE(SPLIT_DIPLOID_GFF.out.hap1, SPLIT_DIPLOID_GFF.out.hap2, fasta)
    }
    else {
        if (!params.source_fasta || !params.source_gff3) {
            error 'Provide --gff3 for an annotated diploid assembly, or --source_fasta and --source_gff3 for Liftoff'
        }
        PREPARE_REFERENCE(fasta)
        def liftoff_input = PREPARE_REFERENCE.out.hap1_fasta
            .combine(PREPARE_REFERENCE.out.hap1_map)
            .map { target, chroms -> tuple('hap1', target, chroms) }
            .mix(PREPARE_REFERENCE.out.hap2_fasta
                .combine(PREPARE_REFERENCE.out.hap2_map)
                .map { target, chroms -> tuple('hap2', target, chroms) })
        LIFTOFF_ANNOTATE(
            liftoff_input,
            requiredFile(params.source_fasta, 'source_fasta'),
            requiredFile(params.source_gff3, 'source_gff3')
        )
        BUILD_ANNOTATION_BUNDLE(
            haplotypeAnnotation(LIFTOFF_ANNOTATE.out.annotation, 'hap1'),
            haplotypeAnnotation(LIFTOFF_ANNOTATE.out.annotation, 'hap2'),
            fasta)
    }

    emit:
    summary        = BUILD_ANNOTATION_BUNDLE.out.summary
    bed            = BUILD_ANNOTATION_BUNDLE.out.bed
    saf            = BUILD_ANNOTATION_BUNDLE.out.saf
    tx2gene        = BUILD_ANNOTATION_BUNDLE.out.tx2gene
    transcriptome  = BUILD_ANNOTATION_BUNDLE.out.transcriptome
    hap1_tx        = BUILD_ANNOTATION_BUNDLE.out.hap1_transcriptome
    hap2_tx        = BUILD_ANNOTATION_BUNDLE.out.hap2_transcriptome
}


workflow DNA_WF {
    take:
    samplesheet
    fasta
    genes_saf
    genes_bed
    annotation_summary

    main:
    def input = samplesheet.splitCsv(header: true, sep: '\t')
        .filter { row -> row.dna }
        .map { row -> tuple(
            [sample: row.dna_id, dna_id: row.dna_id, condition: row.condition],
            row.dna.split(';').collect { path -> file(path, checkIfExists: true) }
        ) }
        .unique { item -> item[0].dna_id }
        .branch { _meta, paths ->
            aligned: paths.size() == 1 && paths[0].name.toLowerCase().endsWith('.bam')
            raw: true
        }
    NORMALIZE_DNA_READS(input.raw)
    DNA_ALIGN(NORMALIZE_DNA_READS.out.reads, fasta)
    PREPARE_DNA_BAM(input.aligned.map { meta, paths -> tuple(meta, paths[0]) })
    def bam = DNA_ALIGN.out.bam.mix(PREPARE_DNA_BAM.out.bam)
    DNA_MAX_SCORE(bam)
    DNA_FEATURE_ASSIGN(DNA_MAX_SCORE.out.maxscore, genes_saf)
    DNA_SPLIT(DNA_FEATURE_ASSIGN.out.split_input)
    DNA_COVERAGE(DNA_SPLIT.out.bams, genes_bed)
    COMPUTE_COPY_NUMBER(
        samplesheet, annotation_summary,
        DNA_COVERAGE.out.coverage.map { _meta, paths -> paths }.flatten().collect(),
        params.min_haplotype_depth, params.min_gene_depth, params.gene_copies,
        params.cn_change_threshold, params.copy_number_level
    )

    emit:
    copy_number = COMPUTE_COPY_NUMBER.out.copy_number
}


workflow RNA_WF {
    take:
    samplesheet
    transcriptome
    tx2gene

    main:
    def input = samplesheet.splitCsv(header: true, sep: '\t')
        .filter { row -> row.rna }
        .map { row -> tuple(
            [sample: row.sample, condition: row.condition, dna_id: row.dna_id],
            row.rna.split(';').collect { path -> file(path, checkIfExists: true) }
        ) }
        .branch { _meta, paths ->
            aligned: paths.size() == 1 && paths[0].name.toLowerCase().endsWith('.bam')
            raw: true
        }
    NORMALIZE_RNA_READS(input.raw)
    RNA_ALIGN(NORMALIZE_RNA_READS.out.reads, transcriptome)
    PREPARE_RNA_BAM(input.aligned.map { meta, paths -> tuple(meta, paths[0]) })
    def bam = RNA_ALIGN.out.bam.mix(PREPARE_RNA_BAM.out.bam)
    OARFISH_ASSIGN(bam, tx2gene, params.oarfish_score, params.oarfish_display_threshold,
                   params.oarfish_strand)
    RNA_HAPLOTYPE_COUNT(OARFISH_ASSIGN.out.assignments,
                        params.gene_resolve_probability, params.isoform_resolve_probability,
                        params.gene_min_probability, params.isoform_min_probability)
    MERGE_RNA_COUNTS(RNA_HAPLOTYPE_COUNT.out.gene_counts.map { _meta, path -> path }.collect(),
                     'gene_counts.tsv')
    MERGE_ISOFORM_COUNTS(RNA_HAPLOTYPE_COUNT.out.isoform_counts.map { _meta, path -> path }.collect(),
                         'isoform_counts.tsv')

    emit:
    counts         = MERGE_RNA_COUNTS.out.counts
    isoform_counts = MERGE_ISOFORM_COUNTS.out.counts
}


workflow PRIORS_WF {
    take:
    transcriptome
    tx2gene
    hap1_tx
    hap2_tx
    isoform_usage

    main:
    def tiled = hap1_tx
        .map { fasta -> tuple([sample: 'SIM_HAP1'], 'hap1', fasta) }
        .mix(hap2_tx.map { fasta -> tuple([sample: 'SIM_HAP2'], 'hap2', fasta) })
    TILE_TRANSCRIPTS(tiled, isoform_usage, tx2gene, params.tiling_read_length,
                     params.tiling_step, params.tiling_max_replicates)
    SIM_RNA_ALIGN(TILE_TRANSCRIPTS.out.reads, transcriptome)
    SIM_OARFISH_ASSIGN(SIM_RNA_ALIGN.out.bam, tx2gene, params.oarfish_score,
                       params.oarfish_display_threshold, params.oarfish_strand)
    SIM_RNA_HAPLOTYPE_COUNT(SIM_OARFISH_ASSIGN.out.assignments,
                            params.gene_resolve_probability, params.isoform_resolve_probability,
                            params.gene_min_probability, params.isoform_min_probability)
    BUILD_MAPPING_PRIORS(
        simulatedCounts(SIM_RNA_HAPLOTYPE_COUNT.out.gene_counts, 'SIM_HAP1'),
        simulatedCounts(SIM_RNA_HAPLOTYPE_COUNT.out.gene_counts, 'SIM_HAP2'),
        simulatedCounts(SIM_RNA_HAPLOTYPE_COUNT.out.isoform_counts, 'SIM_HAP1'),
        simulatedCounts(SIM_RNA_HAPLOTYPE_COUNT.out.isoform_counts, 'SIM_HAP2'),
        TILE_TRANSCRIPTS.out.weighting.first(),
        params.min_simulated_reads)

    emit:
    priors         = BUILD_MAPPING_PRIORS.out.priors
    isoform_priors = BUILD_MAPPING_PRIORS.out.isoform_priors
}


workflow DIFFASE_WF {
    take:
    ase_input

    main:
    COMPILE_ASE_MODEL(channel.value(file("${projectDir}/bin/ASE_model.stan", checkIfExists: true)))
    SPLIT_ASE_INPUT(ase_input, params.contrast, params.ase_shards)
    RUN_ASE_SHARD(SPLIT_ASE_INPUT.out.shards.flatten(), COMPILE_ASE_MODEL.out.model,
                  params.ase_iterations, params.ase_warmup, params.seed)
    MERGE_ASE_RESULTS(RUN_ASE_SHARD.out.results.collect())

    emit:
    results = MERGE_ASE_RESULTS.out.results
}


workflow {
    validateParams()
    log.info "longBayesASE-CN ${workflow.manifest.version ?: 'development'} | step: ${params.step}"

    if (!(params.step in ['all', 'annotation', 'dna', 'rna', 'diffase'])) {
        error '--step must be one of: all, annotation, dna, rna, diffase'
    }
    if (!(params.level in ['gene', 'isoform'])) {
        error "--level must be 'gene' or 'isoform'"
    }
    if (params.step in ['all', 'diffase'] && !params.contrast) {
        error 'Missing required parameter --contrast (example: Treatment:Control, or Case1:Control,Case2:Control)'
    }

    if (params.step == 'annotation') {
        REFERENCE_WF(requiredFile(params.fasta, 'fasta'))
    }
    else if (params.step == 'dna') {
        def sheet = VALIDATE_SAMPLESHEET(requiredFile(params.samplesheet, 'samplesheet'), 'dna').samplesheet
        def fasta = requiredFile(params.fasta, 'fasta')
        REFERENCE_WF(fasta)
        DNA_WF(sheet, fasta, REFERENCE_WF.out.saf, REFERENCE_WF.out.bed, REFERENCE_WF.out.summary)
    }
    else if (params.step == 'rna') {
        def sheet = VALIDATE_SAMPLESHEET(requiredFile(params.samplesheet, 'samplesheet'), 'rna').samplesheet
        REFERENCE_WF(requiredFile(params.fasta, 'fasta'))
        RNA_WF(sheet, REFERENCE_WF.out.transcriptome, REFERENCE_WF.out.tx2gene)
    }
    else if (params.step == 'diffase') {
        if (!params.counts && !params.rna_counts) error 'Provide --counts or --rna_counts'
        def sheet = VALIDATE_SAMPLESHEET(requiredFile(params.samplesheet, 'samplesheet'), 'diffase').samplesheet
        // The transcriptome is only rebuilt when the priors have to be simulated.
        // Everything else this step needs from the annotation is a plain map of
        // genes, transcripts and chromosomes, which comes straight from the GFF3.
        def simulate_priors = !params.priors
        // Copy number is always keyed by gene, so isoform-level counts need the
        // transcript-to-gene map; samplesheet ploidy needs the gene-to-chromosome
        // map. Both come from the GFF3 alone.
        def needs_ploidy_map = !params.counts && !params.copy_number
        def needs_tx2gene = isoformLevel() && !params.counts
        def gene_map = null
        def tx2gene = emptyTx2gene()
        if (simulate_priors) {
            REFERENCE_WF(requiredFile(params.fasta, 'fasta'))
            gene_map = REFERENCE_WF.out.summary
            tx2gene = REFERENCE_WF.out.tx2gene
        }
        else if (needs_ploidy_map || needs_tx2gene || params.gff3) {
            if (!params.gff3) {
                error(needs_tx2gene ?
                    'Provide --gff3: isoform-level counts need the transcript-to-gene map to find their copy number' :
                    'Provide --gff3: samplesheet ploidy needs the gene-to-chromosome map')
            }
            BUILD_FEATURE_MAP(requiredFile(params.gff3, 'gff3'))
            gene_map = BUILD_FEATURE_MAP.out.genes
            tx2gene = BUILD_FEATURE_MAP.out.tx2gene
        }
        def priors
        if (params.priors) priors = requiredFile(params.priors, 'priors')
        else {
            PRIORS_WF(REFERENCE_WF.out.transcriptome, REFERENCE_WF.out.tx2gene,
                      REFERENCE_WF.out.hap1_tx, REFERENCE_WF.out.hap2_tx, isoformUsage())
            priors = isoformLevel() ? PRIORS_WF.out.isoform_priors : PRIORS_WF.out.priors
        }
        if (params.counts) {
            PREPARE_COMBINED_ASE_INPUT(sheet, requiredFile(params.counts, 'counts'), priors)
            DIFFASE_WF(PREPARE_COMBINED_ASE_INPUT.out.ase_input)
        }
        else {
            def copy_number
            if (params.copy_number) copy_number = requiredFile(params.copy_number, 'copy_number')
            else {
                COPY_NUMBER_FROM_PLOIDY(sheet, gene_map)
                copy_number = COPY_NUMBER_FROM_PLOIDY.out.copy_number
            }
            PREPARE_ASE_INPUT(sheet, requiredFile(params.rna_counts, 'rna_counts'),
                              copy_number, priors, tx2gene)
            DIFFASE_WF(PREPARE_ASE_INPUT.out.ase_input)
        }
    }
    else {
        def sheet = VALIDATE_SAMPLESHEET(requiredFile(params.samplesheet, 'samplesheet'), 'all').samplesheet
        def fasta = requiredFile(params.fasta, 'fasta')
        REFERENCE_WF(fasta)
        RNA_WF(sheet, REFERENCE_WF.out.transcriptome, REFERENCE_WF.out.tx2gene)
        def copy_number
        if (params.copy_number) copy_number = requiredFile(params.copy_number, 'copy_number')
        else if (samplesHaveDNA(params.samplesheet)) {
            DNA_WF(sheet, fasta, REFERENCE_WF.out.saf, REFERENCE_WF.out.bed, REFERENCE_WF.out.summary)
            copy_number = DNA_WF.out.copy_number
        }
        else {
            COPY_NUMBER_FROM_PLOIDY(sheet, REFERENCE_WF.out.summary)
            copy_number = COPY_NUMBER_FROM_PLOIDY.out.copy_number
        }
        def priors
        if (params.priors) priors = requiredFile(params.priors, 'priors')
        else {
            // The isoform table is already there, so the priors are weighted by
            // expression without the user having to supply anything.
            def usage = RNA_WF.out.isoform_counts
            if (isoformLevel()) usage = isoformUsage()
            else if (params.isoform_usage) usage = requiredFile(params.isoform_usage, 'isoform_usage')
            PRIORS_WF(REFERENCE_WF.out.transcriptome, REFERENCE_WF.out.tx2gene,
                      REFERENCE_WF.out.hap1_tx, REFERENCE_WF.out.hap2_tx, usage)
            priors = isoformLevel() ? PRIORS_WF.out.isoform_priors : PRIORS_WF.out.priors
        }
        def rna_counts = isoformLevel() ? RNA_WF.out.isoform_counts : RNA_WF.out.counts
        PREPARE_ASE_INPUT(sheet, rna_counts, copy_number, priors, REFERENCE_WF.out.tx2gene)
        DIFFASE_WF(PREPARE_ASE_INPUT.out.ase_input)
    }
}
