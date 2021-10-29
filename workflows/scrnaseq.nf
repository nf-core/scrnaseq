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

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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

workflow SCRNASEQ {
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
