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
    params.star_index, params.txp2gene
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { INPUT_CHECK }       from '../subworkflows/local/input_check'
include { KALLISTO_BUSTOOLS } from '../subworkflows/local/kallisto_bustools'
// include { SCRNASEQ_ALEVIN } from '../subworkflows/local/alevin'
// include { STARSOLO } from '../subworkflows/local/starsolo'


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
(protocol, chemistry) = WorkflowScrnaseq.formatProtocol(params.protocol, params.aligner)

ch_input = file(params.input)
ch_genome_fasta = params.genome_fasta ? file(params.genome_fasta) : []
ch_gtf = params.gtf ? file(params.gtf) : []
ch_kallisto_index = params.kallisto_index ? file(kallisto_index) : []
ch_txp2gene = params.txp2gene ? file(txp2gene) : []
kb_workflow = params.kb_workflow

workflow SCRNASEQ {

    /*
    * Check input files and stage input data
    */
    ch_fastq = INPUT_CHECK( ch_input )
        .reads
        .map {
            meta, reads -> meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, reads ]
        }
        .groupTuple(by: [0])
        .map { it -> [ it[0], it[1].flatten() ] }

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
    }

    // // Run salmon alevin pipeline
    // if (params.aligner == "alevin") {
    //     SCRNASEQ_ALEVIN()
    // }

    // // Run STARSolo pipeline
    // if (params.aligner == "star") {
    //     STARSOLO()
    // }

    // TODO multiqc

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
