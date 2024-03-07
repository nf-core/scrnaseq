include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { FASTQC_CHECK } from '../subworkflows/local/fastqc'
include { KALLISTO_BUSTOOLS } from '../subworkflows/local/kallisto_bustools'
include { SCRNASEQ_ALEVIN } from '../subworkflows/local/alevin'
include { STARSOLO } from '../subworkflows/local/starsolo'
include { CELLRANGER_ALIGN } from "../subworkflows/local/align_cellranger"
include { CELLRANGERARC_ALIGN } from "../subworkflows/local/align_cellrangerarc"
include { UNIVERSC_ALIGN } from "../subworkflows/local/align_universc"
include { MTX_CONVERSION } from "../subworkflows/local/mtx_conversion"
include { GTF_GENE_FILTER } from '../modules/local/gtf_gene_filter'

include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'


workflow SCRNASEQ {

    take:
    ch_fastq

    main:

    protocol_config = Utils.getProtocol(workflow, log, params.aligner, params.protocol)
    if (protocol_config['protocol'] == 'auto' && params.aligner != "cellranger") {
        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
    }

    // general input and params
    ch_genome_fasta = Channel.value(params.fasta ? file(params.fasta) : [])
    ch_gtf = params.gtf ? file(params.gtf) : []
    ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
    ch_motifs = params.motifs ? file(params.motifs) : []
    ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config) : []
    ch_txp2gene = params.txp2gene ? file(params.txp2gene) : []
    ch_multiqc_files = Channel.empty()
    if (params.barcode_whitelist) {
        ch_barcode_whitelist = file(params.barcode_whitelist)
    } else if (protocol_config.containsKey("whitelist")) {
        ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}")
    } else {
        ch_barcode_whitelist = []
    }

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = file(params.input)

    //kallisto params
    ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
    kb_workflow = params.kb_workflow

    //salmon params
    ch_salmon_index = params.salmon_index ? file(params.salmon_index) : []

    //star params
    ch_star_index = params.star_index ? file(params.star_index) : []
    star_feature = params.star_feature

    //cellranger params
    ch_cellranger_index = params.cellranger_index ? file(params.cellranger_index) : []

    //universc params
    ch_universc_index = params.universc_index ? file(params.universc_index) : []

    ch_versions     = Channel.empty()
    ch_mtx_matrices = Channel.empty()

    // Run FastQC
    if (!params.skip_fastqc) {
        FASTQC_CHECK ( ch_fastq )
        ch_versions       = ch_versions.mix(FASTQC_CHECK.out.fastqc_version)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CHECK.out.fastqc_zip.map{ meta, it -> it })
    }

    ch_filter_gtf = GTF_GENE_FILTER ( ch_genome_fasta, ch_gtf ).gtf

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        KALLISTO_BUSTOOLS(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_kallisto_index,
            ch_txp2gene,
            protocol_config['protocol'],
            kb_workflow,
            ch_fastq
        )
        ch_versions = ch_versions.mix(KALLISTO_BUSTOOLS.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(KALLISTO_BUSTOOLS.out.counts)
        ch_txp2gene = KALLISTO_BUSTOOLS.out.txp2gene
    }

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        SCRNASEQ_ALEVIN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_transcript_fasta,
            ch_salmon_index,
            ch_txp2gene,
            ch_barcode_whitelist,
            protocol_config['protocol'],
            ch_fastq
        )
        ch_versions = ch_versions.mix(SCRNASEQ_ALEVIN.out.ch_versions)
        ch_multiqc_files = ch_multiqc_files.mix(SCRNASEQ_ALEVIN.out.alevin_results.map{ meta, it -> it })
        ch_mtx_matrices = ch_mtx_matrices.mix(SCRNASEQ_ALEVIN.out.alevin_results)
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        STARSOLO(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_star_index,
            protocol_config['protocol'],
            ch_barcode_whitelist,
            ch_fastq,
            star_feature,
            protocol_config.get('extra_args', ""),
        )
        ch_versions = ch_versions.mix(STARSOLO.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(STARSOLO.out.star_counts)
        ch_star_index = STARSOLO.out.star_index
        ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.for_multiqc)
    }

    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        CELLRANGER_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_cellranger_index,
            ch_fastq,
            protocol_config['protocol']
        )
        ch_versions = ch_versions.mix(CELLRANGER_ALIGN.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(CELLRANGER_ALIGN.out.cellranger_out)
        ch_star_index = CELLRANGER_ALIGN.out.star_index
        ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER_ALIGN.out.cellranger_out.map{
            meta, outs -> outs.findAll{ it -> it.name == "web_summary.html"}
        })
    }

    // Run universc pipeline
    if (params.aligner == "universc") {
        UNIVERSC_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_universc_index,
            protocol_config['protocol'],
            ch_fastq
        )
        ch_versions = ch_versions.mix(UNIVERSC_ALIGN.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(UNIVERSC_ALIGN.out.universc_out)
    }

    // Run cellranger pipeline
    if (params.aligner == "cellrangerarc") {
        CELLRANGERARC_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_motifs,
            ch_cellranger_index,
            ch_fastq,
            ch_cellrangerarc_config
        )
        ch_versions = ch_versions.mix(CELLRANGERARC_ALIGN.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(CELLRANGERARC_ALIGN.out.cellranger_arc_out)
    }

    // Run mtx to h5ad conversion subworkflow
    MTX_CONVERSION (
        ch_mtx_matrices,
        ch_input,
        ch_txp2gene,
        ch_star_index
    )

    //Add Versions from MTX Conversion workflow too
    ch_versions.mix(MTX_CONVERSION.out.ch_versions)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(softwareVersionsToYAML(ch_versions).collectFile(name: 'versions.yml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
