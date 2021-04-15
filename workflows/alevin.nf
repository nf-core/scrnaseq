////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
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
        .into { gtf_extract_transcriptome; gtf_alevin; gtf_gene_map }
}

//Check if TXP2Gene is provided for Alevin
if (!params.gtf && !params.txp2gene && params.aligner == 'alevin'){
  exit 1, "Must provide either a GTF file ('--gtf') or transcript to gene mapping ('--txp2gene') to align with Alevin"
}

//Setup FastA channels
if( params.fasta ){
    Channel
        .fromPath(params.fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
        .into { genome_fasta_extract_transcriptome }
} else {
  genome_fasta_extract_transcriptome = Channel.empty()
}
//Setup Transcript FastA channels
if( params.transcript_fasta ){
  if( params.aligner == "star" && !params.fasta) {
    exit 1, "Transcriptome-only alignment is not valid with the aligner: ${params.aligner}. Transcriptome-only alignment is only valid with '--aligner alevin'"
  }
    Channel
        .fromPath(params.transcript_fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.transcript_fasta}" }
        .into { transcriptome_fasta_alevin }
} else {
  transcriptome_fasta_alevin = Channel.empty()
}

//Setup channel for salmon index if specified
if (params.salmon_index) {
    Channel
        .fromPath(params.salmon_index)
        .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
        .set { salmon_index_alevin }
}

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */

 if(params.input_paths){
         Channel
             .from(params.input_paths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
             .into { read_files_alevin; read_files_star; read_files_kallisto}
     } else {
         Channel
            .fromFilePairs( params.input )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
            .into { read_files_alevin; read_files_star; read_files_kallisto }
}

if (params.txp2gene){
      Channel
      .fromPath(params.txp2gene)
      .set{ ch_txp2gene } 
}

//Whitelist files for STARsolo and Kallisto
whitelist_folder = "$baseDir/assets/whitelist/"

//Automatically set up proper filepaths to the barcode whitelist files bundled with the pipeline
if (params.type == "10x" && !params.barcode_whitelist){
  barcode_filename = "$whitelist_folder/${params.type}_${params.chemistry}_barcode_whitelist.txt.gz"
  Channel.fromPath(barcode_filename)
         .ifEmpty{ exit 1, "Cannot find ${params.type} barcode whitelist: $barcode_filename" }
         .set{ barcode_whitelist_gzipped }
  Channel.empty().into{ barcode_whitelist_alevinqc }
} else if (params.barcode_whitelist){
  Channel.fromPath(params.barcode_whitelist)
         .ifEmpty{ exit 1, "Cannot find ${params.type} barcode whitelist: $barcode_filename" }
         .into{ barcode_whitelist_alevinqc }
  barcode_whitelist_gzipped = Channel.empty()
}

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def salmon_index_options            = modules['salmon_index']
def gffread_txp2gene_options        = modules['gffread_tx2pgene']
def salmon_alevin_options           = modules['salmon_alevin_options']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GFFREAD as GFFREAD_TRANSCRIPTOME }  from './modules/local/gffread_transcriptome' addParams( options: [:] )
include { SALMON_ALEVIN }                     from './modules/local/salmon_alevin'          addParams( options: salmon_alevin_options )
////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../modules/nf-core/software/gunzip/main'              addParams( options: [:] )
include { GFFREAD as GFFREAD_TXP2GENE } from '../modules/nf-core/software/gffread/main'             addParams( options: gffread_txp2gene_options )
include { SALMON_INDEX }                from '../modules/nf-core/software/salmon/index/main'        addParams( options: salmon_index_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def multiqc_report    = []

workflow SCRNASEQ_ALEVIN {
    ch_software_versions = Channel.empty()

    // unzip barcodes
    if (params.type == "10x" && !params.barcode_whitelist) {
        GUNZIP( barcode_whitelist_gzipped )
    }

    // Preprocessing - Extract transcriptome fasta from genome fasta
    if (params.transcript_fasta) {
        GFFREAD_TRANSCRIPTOME( genome_fasta_extract_transcriptome, gtf_extract_transcriptome)
    }
    
    // build salmon index if not supplied
    if (!params.salmon_index) {
        SALMON_INDEX( GFFREAD_TRANSCRIPTOME.out.transcriptome_extracted)
        salmon_index_alevin = SALMON_INDEX.out.index
    }

    // build the gene map
    if (!params.txp2gene){
        GFFREAD_TXP2GENE( gtf_alevin )
        ch_txp2gene = GFFREAD_TXP2GENE.out.gtf
    }

    // run alignment with salmon alevin
    protocol = "chromium"
    SALMON_ALEVIN ( input_paths, salmon_index_alevin, ch_txp2gene, protocol )
    // TODO run salmon alevin (PR for module open)
      
      // TODO run alevinQC (PR for module open)




}