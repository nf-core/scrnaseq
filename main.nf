#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/scrnaseq
========================================================================================
 nf-core/scrnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/scrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/scrnaseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.checkCondaChannels(log)
}

// Check AWS batch settings
Checks.awsBatch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostName(workflow, params, log)

// Check genome key exists if provided
Checks.genomeExists(params, log)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        include { SCRNASEQ_ALEVIN } from './workflows/alevin'
        SCRNASEQ_ALEVIN()
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        include { STARSOLO } from './workflows/starsolo'
        STARSOLO()
    }

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        include { KALLISTO_BUSTOOLS } from './workflows/kallisto_bustools'
        KALLISTO_BUSTOOLS()
    }
    
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
