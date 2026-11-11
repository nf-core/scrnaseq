#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/scrnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/scrnaseq
    Website: https://nf-co.re/scrnaseq
    Slack  : https://nfcore.slack.com/channels/scrnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCRNASEQ                } from './workflows/scrnaseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_scrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Params cannot be changed if they have been set beforehand
// Thus, manually provided files are not overwritten by the genome attributes

// As discussed in #371 it is desirable for users to be able to provide indices via
// custom igenomes configs in addition to being able to provide them via params directly
params.fasta                = getGenomeAttribute('fasta')
params.gtf                  = getGenomeAttribute('gtf')
params.gff                  = getGenomeAttribute('gff')
params.star_index           = getGenomeAttribute('star')
params.simpleaf_index       = getGenomeAttribute('simpleaf')
params.kallisto_index       = getGenomeAttribute('kallisto')
params.cellranger_index     = getGenomeAttribute('cellranger')
params.txp2gene             = getGenomeAttribute('txp2gene')
params.transcript_fasta     = getGenomeAttribute('transcript_fasta')
params.motifs               = getGenomeAttribute('motifs')
params.cellranger_vdj_index = getGenomeAttribute('cellranger_vdj')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_SCRNASEQ {

    take:
    samplesheet                  // channel: samplesheet read in from --input
    fasta                        // val: path-like string (or null)
    gtf                          // val: path-like string (or null)
    gff                          // val: path-like string (or null)
    star_index                   // val: path-like string (or null)
    simpleaf_index               // val: path-like string (or null)
    kallisto_index               // val: path-like string (or null)
    cellranger_index             // val: path-like string (or null)
    txp2gene                     // val: path-like string (or null)
    transcript_fasta             // val: path-like string (or null)
    motifs                       // val: path-like string (or null)
    cellranger_vdj_index         // val: path-like string (or null)
    multiqc_config               // val: path-like string (or null)
    multiqc_logo                 // val: path-like string (or null)
    multiqc_methods_description  // val: path-like string (or null)
    outdir                       // val: string

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SCRNASEQ (
        samplesheet,
        fasta,
        gtf,
        gff,
        star_index,
        simpleaf_index,
        kallisto_index,
        cellranger_index,
        txp2gene,
        transcript_fasta,
        motifs,
        cellranger_vdj_index,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description,
        outdir,
    )
    emit:
    multiqc_report = SCRNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SCRNASEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.fasta,
        params.gtf,
        params.gff,
        params.star_index,
        params.simpleaf_index,
        params.kallisto_index,
        params.cellranger_index,
        params.txp2gene,
        params.transcript_fasta,
        params.motifs,
        params.cellranger_vdj_index,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_SCRNASEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
