/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { getGenomeAttribute                                } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { FASTQC_CHECK                                      } from '../subworkflows/local/fastqc'
include { KALLISTO_BUSTOOLS                                 } from '../subworkflows/local/kallisto_bustools'
include { SIMPLEAF                                          } from '../subworkflows/local/simpleaf'
include { STARSOLO                                          } from '../subworkflows/local/starsolo'
include { CELLRANGER_ALIGN                                  } from "../subworkflows/local/align_cellranger"
include { CELLRANGER_MULTI_ALIGN                            } from "../subworkflows/local/align_cellrangermulti"
include { CELLRANGERARC_ALIGN                               } from "../subworkflows/local/align_cellrangerarc"
include { MTX_TO_H5AD                                       } from '../modules/local/mtx_to_h5ad'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA } from '../subworkflows/nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { GTF_GENE_FILTER                                   } from '../modules/local/gtf_gene_filter'
include { GUNZIP as GUNZIP_FASTA                            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { H5AD_CONVERSION                                   } from '../subworkflows/local/h5ad_conversion'
include { CONCATENATE_VDJ                                   } from '../modules/local/concatenate_vdj'
include { CONVERT_MUDATA                                    } from '../modules/local/convert_mudata'


workflow SCRNASEQ {

    take:
    ch_fastq

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions      = Channel.empty()
    ch_mtx_matrices  = Channel.empty()

    protocol_config = Utils.getProtocol(workflow, log, params.aligner, params.protocol)
    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
    }

    // general input and params
    ch_genome_fasta         = params.fasta                ? file(params.fasta, checkIfExists: true)    : []
    ch_gtf                  = params.gtf                  ? file(params.gtf, checkIfExists: true)      : []
    ch_transcript_fasta     = params.transcript_fasta     ? file(params.transcript_fasta)              : []
    ch_motifs               = params.motifs               ? file(params.motifs)                        : []
    ch_txp2gene             = params.txp2gene             ? file(params.txp2gene, checkIfExists: true) : []

    if (params.barcode_whitelist) {
        ch_barcode_whitelist = file(params.barcode_whitelist, checkIfExists: true)
    } else if (protocol_config.containsKey("whitelist")) {
        ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}", checkIfExists: true)
    } else {
        ch_barcode_whitelist = []
    }

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = file(params.input)

    //kallisto params
    ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index, checkIfExists: true) : []
    kb_t1c            = params.kb_t1c         ? file(params.kb_t1c, checkIfExists: true) : []
    kb_t2c            = params.kb_t2c         ? file(params.kb_t2c, checkIfExists: true) : []

    //simpleaf params
    ch_simpleaf_index   = params.simpleaf_index ? file(params.simpleaf_index, checkIfExists: true) : []

    //star params
    star_index        = params.star_index ? file(params.star_index, checkIfExists: true) : null
    ch_star_index     = star_index ? Channel.value( [[id: star_index.baseName], star_index] ) : []

    //cellranger params
    ch_cellranger_index = params.cellranger_index ? file(params.cellranger_index, checkIfExists: true) : []

    //cellrangermulti params
    cellranger_vdj_index = params.cellranger_vdj_index      ? file(params.cellranger_vdj_index, checkIfExists: true)      : []
    ch_multi_samplesheet = params.cellranger_multi_barcodes ? file(params.cellranger_multi_barcodes, checkIfExists: true) : []
    empty_file           = file("$projectDir/assets/EMPTY", checkIfExists: true)

    // cellrangerarc params
    ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config)          : []

    // Run FastQC
    if (!params.skip_fastqc) {
        FASTQC_CHECK ( ch_fastq )
        ch_versions      = ch_versions.mix(FASTQC_CHECK.out.fastqc_version)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CHECK.out.fastqc_multiqc.flatten())
    }

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta) {
        if (params.fasta.endsWith('.gz')) {
            ch_genome_fasta    = GUNZIP_FASTA ( [ [:], ch_genome_fasta ] ).gunzip.map { it[1] }
            ch_versions        = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_genome_fasta = Channel.value( ch_genome_fasta )
        }
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], ch_gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value( ch_gtf )
        }
    }

    // filter gtf
    ch_filter_gtf = ch_gtf ? GTF_GENE_FILTER ( ch_genome_fasta, ch_gtf ).gtf : []

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        KALLISTO_BUSTOOLS(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_kallisto_index,
            ch_txp2gene,
            kb_t1c,
            kb_t2c,
            protocol_config['protocol'],
            params.kb_workflow,
            ch_fastq
        )
        ch_versions = ch_versions.mix(KALLISTO_BUSTOOLS.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix( KALLISTO_BUSTOOLS.out.counts_raw, KALLISTO_BUSTOOLS.out.counts_filtered )
        ch_txp2gene = KALLISTO_BUSTOOLS.out.txp2gene
    }

    // Run simpleaf pipeline
    if ( params.aligner == "simpleaf" ) {

        SIMPLEAF(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_transcript_fasta,
            ch_simpleaf_index,
            ch_txp2gene,
            ch_barcode_whitelist,
            protocol_config['protocol'],
            params.simpleaf_umi_resolution,
            ch_fastq,
            [] // for existing map dir; not applicable
        )
        ch_versions = ch_versions.mix(SIMPLEAF.out.ch_versions)
        ch_multiqc_files = ch_multiqc_files.mix(SIMPLEAF.out.quant.map{ _meta, it -> it })
        ch_mtx_matrices = ch_mtx_matrices.mix(
            SIMPLEAF.out.quant.map{
                meta, files -> [
                    meta +
                    [input_type: meta["filtered"] ? "filtered" : "raw" ],
                    files
                ]
            }
        )
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        STARSOLO(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_star_index,
            protocol_config['protocol'],
            ch_barcode_whitelist,
            ch_fastq,
            params.star_feature,
            protocol_config.get('extra_args', ""),
        )
        ch_versions = ch_versions.mix(STARSOLO.out.ch_versions)
        ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.for_multiqc)
        ch_mtx_matrices = ch_mtx_matrices.mix( STARSOLO.out.raw_counts, STARSOLO.out.filtered_counts )
    }

    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        CELLRANGER_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_cellranger_index,
            ch_fastq,
            protocol_config['protocol']
        )
        ch_versions = ch_versions.mix(CELLRANGER_ALIGN.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_ALIGN.out.cellranger_matrices_raw, CELLRANGER_ALIGN.out.cellranger_matrices_filtered )
        ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER_ALIGN.out.cellranger_out.map {
            meta, outs -> outs.findAll{ it -> it.name == "web_summary.html"}
        })
    }

    // Run cellrangerarc pipeline
    if (params.aligner == "cellrangerarc") {
        CELLRANGERARC_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_motifs,
            ch_cellranger_index,
            ch_fastq,
            ch_cellrangerarc_config
        )
        ch_versions = ch_versions.mix(CELLRANGERARC_ALIGN.out.ch_versions)
        ch_mtx_matrices = ch_mtx_matrices.mix(CELLRANGERARC_ALIGN.out.cellranger_arc_out)
    }

    // Run cellrangermulti pipeline
    if (params.aligner == 'cellrangermulti') {

        // parse the input data to generate a collected channel per sample, which will have
        // the metadata and data for each data-type of every sample.
        // then, inside the subworkflow, it can be parsed to manage inputs to the module
        ch_fastq
        .map { meta, fastqs ->
            def parsed_meta = meta.clone() + [ "${meta.feature_type.toString()}": fastqs ]
            parsed_meta.options = [:]

            // add an universal key to differentiate from empty channels so that the "&& meta_gex?.options" lines in the module main.nf can work properly
            parsed_meta.options['data-available'] = true

            // add cellranger options that are currently handled by pipeline, coming from samplesheet
            // the module parses them from the 'gex' options
            if (meta.feature_type.toString() == 'gex') {
                parsed_meta.options['create-bam'] = params.save_align_intermeds  // force bam creation -- param required by cellranger multi
                if (meta.expected_cells) { parsed_meta.options['expected-cells'] = meta.expected_cells }
            }

            [ parsed_meta.id , parsed_meta ]
        }
        .groupTuple( by: 0 )
        .map{ sample_id, map_collection ->
            // Now we must check if every data possibility taken into account in the .branch() operation
            // performed inside the CELLRANGER_MULTI_ALIGN subworkflow are initialized, even with empty files
            // This to ensure that the sizes of each data channel is the same, and the the order and the data types
            // are used together with its rightful pairs
            //
            // data.types: gex, vdj, ab, beam, crispr, cmo

            // clone ArrayBag (received from .groupTuple()) to avoid mutating the input
            def map_collection_clone = []
            map_collection_clone.addAll(map_collection)

            // generate the expected EMPTY tuple when a data type is not used
            // needs to have a collected map like that, so every sample from the samplesheet is analysed one at a time,
            // allowing to have multiple samples in the sheet, having all the data-type tuples initialized,
            // either empty or populated. It will be branched inside the subworkflow.
            if (!map_collection_clone.any{ it.feature_type == 'gex' })    { map_collection_clone.add( [id: sample_id, feature_type: 'gex'   , gex:    empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ it.feature_type == 'vdj' })    { map_collection_clone.add( [id: sample_id, feature_type: 'vdj'   , vdj:    empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ it.feature_type == 'ab' })     { map_collection_clone.add( [id: sample_id, feature_type: 'ab'    , ab:     empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ it.feature_type == 'beam' })   { map_collection_clone.add( [id: sample_id, feature_type: 'beam'  , beam:   empty_file, options:[:] ] ) } // currently not implemented, the input samplesheet checking will not allow it.
            if (!map_collection_clone.any{ it.feature_type == 'crispr' }) { map_collection_clone.add( [id: sample_id, feature_type: 'crispr', crispr: empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ it.feature_type == 'cmo' })    { map_collection_clone.add( [id: sample_id, feature_type: 'cmo'   , cmo:    empty_file, options:[:] ] ) }

            // return final map
            map_collection_clone
        }
        .set{ ch_cellrangermulti_collected_channel }

        // Run cellranger multi
        CELLRANGER_MULTI_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_cellrangermulti_collected_channel,
            ch_cellranger_index,
            cellranger_vdj_index,
            ch_multi_samplesheet
        )
        ch_versions = ch_versions.mix(CELLRANGER_MULTI_ALIGN.out.ch_versions)
        ch_multiqc_files = ch_multiqc_files.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_out.map{
            meta, outs -> outs.findAll{ it -> it.name == "web_summary.html" }
        })
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_raw, CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_filtered )

    }

    //
    // MODULE: Convert mtx matrices to h5ad
    //
    MTX_TO_H5AD (
        ch_mtx_matrices,
        ch_txp2gene,
        star_index ? ch_star_index.map{it[1]} : [],
        params.aligner
    )
    ch_versions = ch_versions.mix(MTX_TO_H5AD.out.versions.first())
    ch_h5ads = MTX_TO_H5AD.out.h5ad

    //
    // SUBWORKFLOW: Run cellbender remove background subworkflow
    //
    if ( !params.skip_cellbender && !(params.aligner in ['cellrangerarc']) ) {
        // module should only run on the raw matrices thus, filter-out the filtered result of the aligners that can produce it
        H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA (
            ch_h5ads
                .filter { meta, mtx_files -> meta.input_type == 'raw' }
                .map { meta, mtx_files -> [ meta + [input_type: 'cellbender_filter'], mtx_files ]} // to avoid name collision
        )
        ch_h5ads = ch_h5ads.mix(
            H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA.out.h5ad
        )
    }

    //
    // SUBWORKFLOW: Concat samples and convert h5ad to other formats
    //
    H5AD_CONVERSION (
        ch_h5ads,
        ch_input
    )
    ch_versions = ch_versions.mix(H5AD_CONVERSION.out.ch_versions)
    //
    // MODULE: Concat vdj samples and save as h5ad format
    //
    if (params.aligner == "cellrangermulti") {
        CONCATENATE_VDJ (
            CELLRANGER_MULTI_ALIGN.out.vdj
        )
        ch_versions = ch_versions.mix(CONCATENATE_VDJ.out.versions)

    //
    // SUBWORKFLOW: Concat GEX, VDJ and CITE data and save as MuData object
    //

        ch_vdj = CONCATENATE_VDJ.out.h5ad
            .map { meta, file -> [meta, file] }
            .ifEmpty { [[id: 'dummy'], []] }
    } else {
        ch_vdj = [[id: 'dummy'], []]
    }

    ch_h5ad_concat = H5AD_CONVERSION.out.h5ads

    // Filter input_type:'filtered'
    ch_h5ad_concat_filtered = ch_h5ad_concat.filter { item ->
        item[0].input_type == 'filtered'
    }

    if (params.aligner == "cellrangermulti") {
        CONVERT_MUDATA(
            ch_h5ad_concat_filtered,
            ch_vdj
        )
        ch_versions = ch_versions.mix(CONVERT_MUDATA.out.versions)
    } else {'nothing to convert to MuData'}


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'scrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
