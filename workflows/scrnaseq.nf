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
include { gtfSourceFixNeeded                                } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { isStarIndexLegacy                                 } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
include { PREPARE_GENOME                                    } from '../subworkflows/local/prepare_genome'
include { FASTQC_CHECK                                      } from '../subworkflows/local/fastqc'
include { KALLISTO_BUSTOOLS                                 } from '../subworkflows/local/kallisto_bustools'
include { SIMPLEAF                                          } from '../subworkflows/local/simpleaf'
include { STARSOLO                                          } from '../subworkflows/local/starsolo'
include { CELLRANGER_ALIGN                                  } from "../subworkflows/local/align_cellranger"
include { CELLRANGER_MULTI_ALIGN                            } from "../subworkflows/local/align_cellrangermulti"
include { CELLRANGERARC_ALIGN                               } from "../subworkflows/local/align_cellrangerarc"
include { MTX_TO_H5AD                                       } from '../modules/local/mtx_to_h5ad'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA } from '../subworkflows/nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { H5AD_CONVERSION                                   } from '../subworkflows/local/h5ad_conversion'


workflow SCRNASEQ {

    take:
    ch_fastq                    // channel: [ meta, fastq ] from samplesheet
    fasta                       // val: path-like string (or null)
    gtf                         // val: path-like string (or null)
    gff                         // val: path-like string (or null)
    star_index                  // val: path-like string (or null)
    simpleaf_index              // val: path-like string (or null)
    kallisto_index              // val: path-like string (or null)
    cellranger_index            // val: path-like string (or null)
    txp2gene                    // val: path-like string (or null)
    transcript_fasta            // val: path-like string (or null)
    motifs                      // val: path-like string (or null)
    cellranger_vdj_index        // val: path-like string (or null)
    multiqc_config              // val: path-like string (or null)
    multiqc_logo                // val: path-like string (or null)
    multiqc_methods_description // val: path-like string (or null)
    outdir                      // val: string

    main:
    ch_multiqc_files = channel.empty()
    ch_versions      = channel.empty()
    ch_mtx_matrices  = channel.empty()

    protocol_config = Utils.getProtocol(workflow, log, params.aligner, params.protocol)
    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
    }

    // Get qcatch chemistry for simpleaf QC (derived from the simpleaf protocol; null if unmapped)
    qcatch_chemistry = params.aligner == "simpleaf" ? protocol_config['qcatch_protocol'] : null

    // general input and params
    ch_transcript_fasta     = transcript_fasta ? file(transcript_fasta, checkIfExists: true) : []
    ch_motifs               = motifs           ? file(motifs, checkIfExists: true)           : []
    ch_txp2gene             = txp2gene         ? file(txp2gene, checkIfExists: true)         : []

    if (params.barcode_whitelist) {
        ch_barcode_whitelist = file(params.barcode_whitelist, checkIfExists: true)
    } else if (protocol_config.containsKey("whitelist")) {
        ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}", checkIfExists: true)
    } else {
        ch_barcode_whitelist = []
    }

    // Warn if both GTF and GFF files are provided
    if (gtf && gff) {
        log.warn("Both GTF and GFF files are provided. GTF file will be used.")
    }

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the
    // AnnData objects.
    ch_input = file(params.input)

    //kallisto params
    ch_kallisto_index = kallisto_index ? file(kallisto_index, checkIfExists: true) : []
    kb_t1c            = params.kb_t1c  ? file(params.kb_t1c, checkIfExists: true)  : []
    kb_t2c            = params.kb_t2c  ? file(params.kb_t2c, checkIfExists: true)  : []

    //simpleaf params
    ch_simpleaf_index   = simpleaf_index ? file(simpleaf_index, checkIfExists: true) : []

    //star params
    star_index_legacy = isStarIndexLegacy(params.genome, params.genomes, star_index) ?: false
    star_index        = star_index ? file(star_index, checkIfExists: true) : null

    //cellranger params
    ch_cellranger_index = cellranger_index ? file(cellranger_index, checkIfExists: true) : []

    //cellrangermulti params
    cellranger_vdj_index = cellranger_vdj_index             ? file(cellranger_vdj_index, checkIfExists: true)             : []
    ch_multi_samplesheet = params.cellranger_multi_barcodes ? file(params.cellranger_multi_barcodes, checkIfExists: true) : []
    empty_file           = file("$projectDir/assets/EMPTY", checkIfExists: true)

    // cellrangerarc params
    ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config)          : []

    // Run FastQC
    if (!params.skip_fastqc) {
        FASTQC_CHECK ( ch_fastq )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CHECK.out.fastqc_multiqc.flatten())
    }

    //
    // Prepare reference FASTA and GTF (gunzip, filter, optional Cell Ranger GTF source fix)
    //
    PREPARE_GENOME(
        fasta,
        gtf,
        gff,
        gtfSourceFixNeeded(params.aligner, params.genome, params.genomes, gtf)
    )
    ch_genome_fasta = PREPARE_GENOME.out.fasta
    ch_filter_gtf   = PREPARE_GENOME.out.gtf
    ch_versions     = ch_versions.mix(PREPARE_GENOME.out.versions)

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
        ch_mtx_matrices = ch_mtx_matrices.mix( KALLISTO_BUSTOOLS.out.counts_raw, KALLISTO_BUSTOOLS.out.counts_filtered )
        ch_txp2gene = KALLISTO_BUSTOOLS.out.txp2gene
        ch_versions = ch_versions.mix(KALLISTO_BUSTOOLS.out.ch_versions)
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
            qcatch_chemistry,
            params.skip_qcatch,
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
            star_index,
            star_index_legacy,
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
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_ALIGN.out.cellranger_matrices_raw, CELLRANGER_ALIGN.out.cellranger_matrices_filtered )
        ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER_ALIGN.out.cellranger_out.map {
            _meta, outs -> outs.findAll{ summary -> summary.name == "web_summary.html"}
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
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGERARC_ALIGN.out.cellrangerarc_mtx_raw, CELLRANGERARC_ALIGN.out.cellrangerarc_mtx_filtered )
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
                parsed_meta.options['chemistry'] = protocol_config['protocol']
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
            if (!map_collection_clone.any{ m -> m.feature_type == 'gex' })    { map_collection_clone.add( [id: sample_id, feature_type: 'gex'   , gex:    empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ m -> m.feature_type == 'vdj' })    { map_collection_clone.add( [id: sample_id, feature_type: 'vdj'   , vdj:    empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ m -> m.feature_type == 'ab' })     { map_collection_clone.add( [id: sample_id, feature_type: 'ab'    , ab:     empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ m -> m.feature_type == 'beam' })   { map_collection_clone.add( [id: sample_id, feature_type: 'beam'  , beam:   empty_file, options:[:] ] ) } // currently not implemented, the input samplesheet checking will not allow it.
            if (!map_collection_clone.any{ m -> m.feature_type == 'crispr' }) { map_collection_clone.add( [id: sample_id, feature_type: 'crispr', crispr: empty_file, options:[:] ] ) }
            if (!map_collection_clone.any{ m -> m.feature_type == 'cmo' })    { map_collection_clone.add( [id: sample_id, feature_type: 'cmo'   , cmo:    empty_file, options:[:] ] ) }

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
        ch_multiqc_files = ch_multiqc_files.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_out.map{
            _meta, outs -> outs.findAll{ it -> it.name == "web_summary.html" }
        })
        ch_mtx_matrices = ch_mtx_matrices.mix( CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_raw, CELLRANGER_MULTI_ALIGN.out.cellrangermulti_mtx_filtered )

    }

    //
    // MODULE: Convert mtx matrices to h5ad
    //
    MTX_TO_H5AD (
        ch_mtx_matrices,
        ch_txp2gene,
        star_index ?: [],
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
                .filter { meta, _mtx_files -> meta.input_type == 'raw' }
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

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'scrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    if (!params.skip_multiqc) {
        //
        // MODULE: MultiQC
        //
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'scrnaseq'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    } else {
        ch_multiqc_report = channel.empty()
    }

    emit:
    multiqc_report = ch_multiqc_report           // channel: [ path(multiqc_report.html) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
