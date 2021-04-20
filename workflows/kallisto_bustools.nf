////////////////////////////////////////////////////
/* --         KALLISTO BUSTOOLS WORKFLOW       -- */
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
        .set { gtf }
}

//Setup FastA channels
if( params.genome_fasta ){
    Channel
        .fromPath(params.genome_fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}" }
        .set { genome_fasta }
} 

//Setup Transcript FastA channels
if( params.transcript_fasta ){
  Channel
        .fromPath(params.transcript_fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.transcript_fasta}" }
        .set { transcriptome_fasta }
} 

// Check if files for index building are given if no index is specified
if (!params.kallisto_index && !params.genome_fasta || !params.transcript_fasta) {
  exit 1, "Must provide a genome fasta file ('--genome_fasta') or a transcript fasta ('--transcript_fasta') if no index is given!"
}

//Setup channel for salmon index if specified
if (params.kallisto_index) {
    Channel
        .fromPath(params.kallisto_index)
        .ifEmpty { exit 1, "Kallisto index not found: ${params.kallisto_index}" }
        .set { ch_kallisto_index }
}

// Kallist gene map
// Check if txp2gene file has been provided
if (params.kallisto_gene_map){
      Channel
      .fromPath(params.kallisto_gene_map)
      .set{ ch_kallisto_gene_map } 
}
if (!params.gtf && !params.kallisto_gene_map){
  exit 1, "Must provide either a GTF file ('--gtf') or kallisto gene mao('--kallisto_gene_map') to align with kallisto bustools!"
}

// Create a channel for input read files
if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

//Whitelist files for STARsolo and Kallisto
whitelist_folder = "$baseDir/assets/whitelist/"

// Get the protocol parameter
protocol = params.protocol

//Automatically set up proper filepaths to the barcode whitelist files bundled with the pipeline
if ((params.protocol == "chromium" || params.protocol == "chromiumV3") && !params.barcode_whitelist){
  barcode_filename = "$whitelist_folder/10x_${params.chemistry}_barcode_whitelist.txt.gz"
  Channel.fromPath(barcode_filename)
         .ifEmpty{ exit 1, "Cannot find ${params.protocol} barcode whitelist: $barcode_filename" }
         .set{ barcode_whitelist_gzipped }
} else if (params.barcode_whitelist){
  Channel.fromPath(params.barcode_whitelist)
         .ifEmpty{ exit 1, "Cannot find ${params.protocol} barcode whitelist: $barcode_filename" }
         .set{ ch_barcode_whitelist }
}

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def kallisto_index_options          = modules['kallistobustools_ref']
def kallistobustools_count_options  = modules['kallistobustools_count']
def multiqc_options                 = modules['multiqc_kb']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { INPUT_CHECK        }                from '../subworkflows/local/input_check'        addParams( options: [:] )
include { GENE_MAP }                          from '../modules/local/gene_map'                addParams( options: [:] )
include { KALLISTOBUSTOOLS_COUNT }            FROM '../modules/local/kallistobustools_count'  addParams( options: kallistobustools_count_options )
include { GET_SOFTWARE_VERSIONS }             from '../modules/local/get_software_versions'   addParams( options: [publish_files: ['csv':'']]       )
include { MULTIQC }                           from '../modules/local/multiqc_kb'              addParams( options: multiqc_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../modules/nf-core/software/gunzip/main'                    addParams( options: [:] )
include { KALLISTOBUSTOOLS_REF }       from '../modules/nf-core/software/kallistobustools/ref/main'      addParams( options: kallistobustools_ref_options )
////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def multiqc_report    = []

workflow KALLISTO_BUSTOOLS {
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
    if ((params.protocol == "chromium" || params.protocol == "chromiumV3") && !params.barcode_whitelist) {
        GUNZIP( barcode_whitelist_gzipped )
        ch_barcode_whitelist = GUNZIP.out.gunzip
    }

    /*
    * Generate Kallisto Gene Map if not supplied
    */ 
    if (!params.kallisto_gene_map) {
      GENE_MAP)( gtf )
      ch_kallisto_gene_map = GENE_MAP.out.gene_map
    }

    /*
    * Generate kallisto index
    */ 
    if (!params.kallisto_index) { 
      val kb_workflow = "standard"
      KALLISTOBUSTOOLS_REF( genome_fasta, gtf, kb_workflow )
    }


    // collect software versions
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect() )

    /*
    * MultiQC
    */
    // if (!params.skip_multiqc) {
    //     workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //         ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }

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