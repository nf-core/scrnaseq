/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScrnaseq.initialise(params, log)

def checkPathParamList = [
    params.input, params.multiqc_config, params.genome_fasta, params.gtf,
    params.transcript_fasta, params.salmon_index, params.kallisto_index,
    params.star_index, params.txp2gene, params.barcode_whitelist
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { KALLISTO_BUSTOOLS } from '../subworkflows/local/kallisto_bustools'
include { SCRNASEQ_ALEVIN } from '../subworkflows/local/alevin'
include { STARSOLO } from '../subworkflows/local/starsolo'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC } from "../modules/local/multiqc"

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
(protocol, chemistry) = WorkflowScrnaseq.formatProtocol(params.protocol, params.aligner)

ch_input = file(params.input)
ch_genome_fasta = params.genome_fasta ? file(params.genome_fasta) : []
ch_gtf = params.gtf ? file(params.gtf) : []
ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
ch_salmon_index = params.salmon_index ? file(params.salmon_index) : []
ch_txp2gene = params.txp2gene ? file(txp2gene) : []
ch_star_index = params.star_index ? file(params.star_index) : []
ch_multiqc_alevin = []
ch_multiqc_star = []
kb_workflow = params.kb_workflow

if (params.barcode_whitelist) {
    ch_barcode_whitelist = file(params.barcode_whitelist)
} else if (params.protocol.contains("10X")) {
    ch_barcode_whitelist = file("$baseDir/assets/whitelist/10x_${chemistry}_barcode_whitelist.txt.gz", checkIfExists: true)
} else {
    ch_barcode_whitelist = []
}

workflow SCRNASEQ {

    ch_versions = Channel.empty()

    // Check input files and stage input data
    ch_fastq = INPUT_CHECK( ch_input )
        .reads
        .map {
            meta, reads -> meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, reads ]
        }
        .groupTuple(by: [0])
        .map { it -> [ it[0], it[1].flatten() ] }

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
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        STARSOLO(
            ch_genome_fasta,
            ch_gtf,
            ch_star_index,
            protocol,
            ch_barcode_whitelist,
            ch_fastq
        )
        ch_versions = ch_versions.mix(STARSOLO.out.ch_versions)
        ch_multiqc_star = STARSOLO.out.for_multiqc
    }

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
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
