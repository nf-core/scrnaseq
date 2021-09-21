////////////////////////////////////////////////////
/* --         STARSolo WORKFLOW                -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

//Check if GTF is supplied properly
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .collect()
        .set { gtf }
}

//Setup FastA channels
if( params.genome_fasta ){
    Channel
        .fromPath(params.genome_fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}" }
        .collect()
        .set { genome_fasta }
} 

//Setup Transcript FastA channels
if( params.transcript_fasta ){
  Channel
        .fromPath(params.transcript_fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.transcript_fasta}" }
        .collect()
        .set { transcriptome_fasta }
} 

// Check if STAR index is supplied properly
if( params.star_index  ){
    star_index = Channel
        .fromPath(params.star_index)
        .collect()
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

if (!params.star_index && (!params.gtf || !params.genome_fasta)){
  exit 1, "STAR needs either a GTF + FASTA or a precomputed index supplied."
}

// Create a channel for input read files
if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }

// Check if txp2gene file has been provided
if (params.txp2gene){
      Channel
      .fromPath(params.txp2gene)
      .collect()
      .set{ ch_txp2gene } 
}

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

// Get the protocol parameter
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, "star")

//Whitelist files for STARsolo and Kallisto
whitelist_folder = "$baseDir/assets/whitelist/"

//Automatically set up proper filepaths to the barcode whitelist files bundled with the pipeline
if (params.protocol.contains("10X") && !params.barcode_whitelist){
  barcode_filename = "$whitelist_folder/10x_${chemistry}_barcode_whitelist.txt.gz"
  Channel.fromPath(barcode_filename)
         .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
         .collect()
         .set{ barcode_whitelist_gzipped }
} else if (params.barcode_whitelist){
  Channel.fromPath(params.barcode_whitelist)
         .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
         .collect()
         .set{ ch_barcode_whitelist }
}

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def star_genomegenerate_options     = modules['star_genomegenerate']
def star_align_options              = modules['star_align']
def multiqc_options                 = modules['multiqc_alevin']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { INPUT_CHECK        }          from '../subworkflows/local/input_check'                    addParams( options: [:] )
include { GET_SOFTWARE_VERSIONS }       from '../modules/local/get_software_versions'               addParams( options: [publish_files: ['csv':'']]       )
include { MULTIQC }                     from '../modules/local/multiqc_alevin'                      addParams( options: multiqc_options )
include { STAR_ALIGN }                  from '../modules/local/star_align'                          addParams( options: star_align_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../modules/nf-core/software/gunzip/main'              addParams( options: [:] )
include { STAR_GENOMEGENERATE }         from '../modules/nf-core/software/star/genomegenerate/main' addParams( options: star_genomegenerate_options )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def multiqc_report    = []

workflow STARSOLO {
    ch_software_versions = Channel.empty()

    /*
    * Check input files and stage input data
    */
    INPUT_CHECK( ch_input )
    .map {
        meta, reads -> meta.id = meta.id.split('_')[0..-2].join('_')
        [ meta, reads ]
    }
    .groupTuple(by: [0])
    .map { it -> [ it[0], it[1].flatten() ] }
    .set { ch_fastq }

    // unzip barcodes
    if (params.protocol.contains("10X") && !params.barcode_whitelist) {
        GUNZIP( barcode_whitelist_gzipped )
        ch_barcode_whitelist = GUNZIP.out.gunzip
    }

    /*
    * Build STAR index if not supplied
    */
    if (!params.star_index) {
      STAR_GENOMEGENERATE( genome_fasta, gtf )
      star_index = STAR_GENOMEGENERATE.out.index
    }

    /*
    * Perform mapping with STAR
    */ 
    STAR_ALIGN( 
      ch_fastq,
      star_index,
      gtf,
      ch_barcode_whitelist,
      protocol
    )
    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.first().ifEmpty(null))
    ch_star_multiqc      = STAR_ALIGN.out.log_final

    // collect software versions
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect() )

    /*
    * MultiQC
    */
    if (!params.skip_multiqc) {
        workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////