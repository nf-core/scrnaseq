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


// Check if files for index building are given if no index is specified
if (!params.kallisto_index && (!params.genome_fasta || !params.gtf)) {
    exit 1, "Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf') if no index is given!"
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
    exit 1, "Must provide either a GTF file ('--gtf') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"
}

// Get the protocol parameter
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, "kallisto")
kb_workflow = "standard"

// Create a channel for input read files
if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def kallistobustools_ref_options    = modules['kallistobustools_ref']
def kallistobustools_count_options  = modules['kallistobustools_count']
def multiqc_options                 = modules['multiqc_kb']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { INPUT_CHECK        }                from '../subworkflows/local/input_check'        addParams( options: [:] )
include { GENE_MAP }                          from '../modules/local/gene_map'                addParams( options: [:] )
include { GET_SOFTWARE_VERSIONS }             from '../modules/local/get_software_versions'   addParams( options: [publish_files: ['csv':'']]       )
include { MULTIQC }                           from '../modules/local/multiqc_kb'              addParams( options: multiqc_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                     from '../modules/nf-core/modules/gunzip/main'                    addParams( options: [:] )
include { KALLISTOBUSTOOLS_REF }       from '../modules/nf-core/modules/kallistobustools/ref/main'       addParams( options: kallistobustools_ref_options )
include { KALLISTOBUSTOOLS_COUNT }     from '../modules/nf-core/modules/kallistobustools/count/main'   addParams( options: kallistobustools_count_options )

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

    /*
    * Generate Kallisto Gene Map if not supplied and index is given
    * If index is given, the gene map will be generated in the 'kb ref' step
    */
    if (!params.kallisto_gene_map && params.kallisto_index) {
        GENE_MAP( gtf )
        ch_kallisto_gene_map = GENE_MAP.out.gene_map
    }

    /*
    * Generate kallisto index
    */
    if (!params.kallisto_index) {
        KALLISTOBUSTOOLS_REF( genome_fasta, gtf, kb_workflow )
        ch_kallisto_gene_map = KALLISTOBUSTOOLS_REF.out.t2g
        ch_kallisto_index    = KALLISTOBUSTOOLS_REF.out.index
    }

    /*
    * Quantification with kallistobustools count
    */
    KALLISTOBUSTOOLS_COUNT(
        ch_fastq,
        ch_kallisto_index.collect(),
        ch_kallisto_gene_map.collect(),
        false,
        false,
        kb_workflow,
        protocol
    )
    ch_software_versions = ch_software_versions.mix(KALLISTOBUSTOOLS_COUNT.out.versions.first().ifEmpty(null))

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
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
