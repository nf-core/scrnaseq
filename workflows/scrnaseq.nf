/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScrnaseq.initialise(params, log)

def checkPathParamList = [
    params.input, params.multiqc_config, params.genome_fasta, params.gtf,
    params.transcript_fasta, params.salmon_index, params.kallisto_index,
    params.star_index, params.txp2gene, params.barcode_whitelist, params.cellranger_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { KALLISTO_BUSTOOLS } from '../subworkflows/local/kallisto_bustools'
include { SCRNASEQ_ALEVIN   } from '../subworkflows/local/alevin'
include { STARSOLO          } from '../subworkflows/local/starsolo'
include { CELLRANGER_ALIGN  } from "../subworkflows/local/align_cellranger"
include { H5AD_CONVERSION   } from "../subworkflows/local/conversion_to_h5ad"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC } from "../modules/local/multiqc"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
// TODO: Are this channels still necessary?
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
(protocol, chemistry, other_parameters) = WorkflowScrnaseq.formatProtocol(params.protocol, params.aligner)

// general input and params
ch_input = file(params.input)
ch_genome_fasta = params.genome_fasta ? file(params.genome_fasta) : []
ch_gtf = params.gtf ? file(params.gtf) : []
ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
ch_txp2gene = params.txp2gene ? file(params.txp2gene) : []
ch_multiqc_alevin = []
ch_multiqc_star = []
if (params.barcode_whitelist) {
    ch_barcode_whitelist = file(params.barcode_whitelist)
} else if (params.protocol.contains("10X")) {
    ch_barcode_whitelist = file("$baseDir/assets/whitelist/10x_${chemistry}_barcode_whitelist.txt.gz", checkIfExists: true)
} else {
    ch_barcode_whitelist = []
}


//kallisto params
ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
kb_workflow = params.kb_workflow

//salmon params
ch_salmon_index = params.salmon_index ? file(params.salmon_index) : []

//star params
ch_star_index = params.star_index ? file(params.star_index) : []

//cellranger params
ch_cellranger_index = params.cellranger_index ? file(params.cellranger_index) : []


workflow SCRNASEQ {

    ch_versions     = Channel.empty()
    ch_mtx_matrices = Channel.empty()

    // Check input files and stage input data
    ch_fastq = INPUT_CHECK( ch_input ).reads

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        KALLISTO_BUSTOOLS(
            ch_genome_fasta,
            ch_gtf,
            ch_kallisto_index,
            ch_txp2gene,
            protocol,
            chemistry,
            kb_workflow,
            ch_fastq
        )
        ch_versions = ch_versions.mix(KALLISTO_BUSTOOLS.out.ch_versions)
    }

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        SCRNASEQ_ALEVIN(
            ch_genome_fasta,
            ch_gtf,
            ch_transcript_fasta,
            ch_salmon_index,
            ch_txp2gene,
            ch_barcode_whitelist,
            protocol,
            chemistry,
            ch_fastq
        )
        ch_versions = ch_versions.mix(SCRNASEQ_ALEVIN.out.ch_versions)
        ch_multiqc_alevin = SCRNASEQ_ALEVIN.out.for_multiqc
        ch_mtx_matrices = ch_mtx_matrices.mix(SCRNASEQ_ALEVIN.out.alevin_results)
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        STARSOLO(
            ch_genome_fasta,
            ch_gtf,
            ch_star_index,
            protocol,
            ch_barcode_whitelist,
            ch_fastq,
            other_parameters
        )
        ch_versions = ch_versions.mix(STARSOLO.out.ch_versions)
        ch_multiqc_star = STARSOLO.out.for_multiqc
    }

    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        CELLRANGER_ALIGN(
            ch_genome_fasta,
            ch_gtf,
            ch_cellranger_index,
            ch_fastq
        )
        ch_versions = ch_versions.mix(CELLRANGER_ALIGN.out.ch_versions)
    }

    // Run mtx to h5ad conversion subworkflow
    H5AD_CONVERSION (
        ch_mtx_matrices,
        ch_input
    )

    // collect software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (!params.skip_multiqc) {
        ch_workflow_summary = Channel.value(
            WorkflowScrnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        MULTIQC(
            ch_multiqc_config,
            ch_multiqc_custom_config,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary,
            ch_multiqc_alevin,
            ch_multiqc_star
        )
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
