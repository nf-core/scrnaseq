/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScrnaseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { SCRNASEQ_ALEVIN } from './alevin'
include { STARSOLO } from './starsolo'
include { KALLISTO_BUSTOOLS } from './kallisto_bustools'


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
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, params.aligner)


workflow SCRNASEQ {

    /*
    * Check input files and stage input data
    */
    INPUT_CHECK( ch_input )
    .reads
    .map {
        meta, reads -> meta.id = meta.id.split('_')[0..-2].join('_')
        [ meta, reads ]
    }
    .groupTuple(by: [0])
    .map { it -> [ it[0], it[1].flatten() ] }
    .set { ch_fastq }

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        SCRNASEQ_ALEVIN()
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        STARSOLO()
    }

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        KALLISTO_BUSTOOLS()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

// Added to individual workflows

/*
========================================================================================
    THE END
========================================================================================
*/
