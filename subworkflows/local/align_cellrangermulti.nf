//
// Include modules
//
include { CELLRANGER_MKGTF                  } from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF                  } from "../../modules/nf-core/cellranger/mkref/main.nf"
include { CELLRANGER_MKVDJREF               } from "../../modules/nf-core/cellranger/mkvdjref/main.nf"
include { CELLRANGER_MULTI                  } from "../../modules/nf-core/cellranger/multi/main.nf"
include { PARSE_CELLRANGERMULTI_SAMPLESHEET } from "../../modules/local/parse_cellrangermulti_samplesheet.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI_ALIGN {
    take:
        ch_fasta
        ch_gtf
        ch_fastq
        cellranger_gex_index
        cellranger_vdj_index
        empty_file
        ch_multi_samplesheet

    main:
        ch_versions    = Channel.empty()

        // since we merged all data as a meta, now we have a channel per sample, which
        // every item is a meta map for each data-type
        // now we can split it back for passing as input to the module
        ch_fastq
        .flatten()
        .map{ meta ->
            def meta_clone = meta.clone()
            def data_dict  = meta_clone.find{ it.key == "${meta_clone.feature_type}" }
            fastqs = data_dict?.value
            meta_clone.remove( data_dict?.key )
            [ meta_clone, fastqs ]
        }
        .branch {
            meta, fastq ->
                gex: meta.feature_type == "gex"
                    return [ meta, fastq ]
                vdj: meta.feature_type == "vdj"
                    return [ meta, fastq ]
                ab: meta.feature_type == "ab"
                    return [ meta, fastq ]
                beam: meta.feature_type == "beam"
                    return [ meta, fastq ]
                crispr: meta.feature_type == "crispr"
                    return [ meta, fastq ]
                cmo: meta.feature_type == "cmo"
                    return [ meta, fastq ]
        }
        .set { ch_grouped_fastq }

        // Assign other cellranger reference files
        ch_gex_frna_probeset      = params.gex_frna_probe_set ? file(params.gex_frna_probe_set) : empty_file
        ch_gex_target_panel       = params.gex_target_panel   ? file(params.gex_target_panel)   : empty_file
        ch_gex_cmo_set            = params.gex_cmo_set        ? file(params.gex_cmo_set)        : empty_file
        ch_gex_barcodes           = params.gex_barcode_sample_assignment ? file(params.gex_barcode_sample_assignment) : empty_file
        ch_fb_reference           = params.fb_reference       ? file(params.fb_reference)       : empty_file
        ch_vdj_primer_index       = params.vdj_inner_enrichment_primers ? file(params.vdj_inner_enrichment_primers) : empty_file
        ch_beam_antigen_panel_csv = empty_file // currently not implemented
        ch_beam_control_panel_csv = empty_file // currently not implemented

        // parse frna and barcode information
        if (ch_multi_samplesheet) {

            //
            // Here, we parse the received cellranger multi barcodes samplesheet.
            // We first use the get the PARSE_CELLRANGERMULTI_SAMPLESHEET module to check it and guarantee structure
            // and also split it to have one fnra/cmo .csv for each sample.
            //
            // The selection of the GEX fastqs is because samples are always expected to have at least GEX data.
            // Then, using "combined" map, which means, the "additional barcode information" of each sample, we then,
            // parse it to generate the cmo / frna samplesheets to be used by each sample.
            //
            // Here, to guarantee it and take advantage of the "FIFO"-rule and are sure that the data used in the
            // module is from the same sample from the "normal" samplesheet. We have to use the .concat().groupTuple()
            // pipe instead of .join() because .join() outputs first the arrays that could be joined and afterwards
            // the ones with "remainders", thus, we would not ensure "FIFO" and the same order.
            //
            // To guarantee this, we can define two nf-tests, one having only one sample with CMO and another with two
            // samples using CMOs, even if wrongly/repeated, but just to guarantee FIFO is working.
            //

            PARSE_CELLRANGERMULTI_SAMPLESHEET( ch_multi_samplesheet )

            ch_grouped_fastq.gex
            .map{ [it[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.cmo.flatten().map { [ "${it.baseName}" - "_cmo", it ] } )
            .groupTuple()
            .map { if ( it.size() == 2 ) { it[1] } else { empty_file } } // a correct tuple from snippet will have: [ sample, cmo.csv ]
            .flatten()
            .set { ch_cmo_barcode_csv }

            ch_grouped_fastq.gex
            .map{ [it[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.frna.flatten().map { [ "${it.baseName}" - "_frna", it ] } )
            .groupTuple()
            .map { if ( it.size() == 2 ) { it[1] } else { empty_file } } // a correct tuple from snippet will have: [ sample, frna.csv ]
            .flatten()
            .set { ch_frna_sample_csv }

        } else {
            ch_cmo_barcode_csv = empty_file
            ch_frna_sample_csv = empty_file
        }

        //
        // Prepare GTF
        //
        if (!cellranger_gex_index || !cellranger_vdj_index) {

            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF ( ch_gtf )
            ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

        }

        //
        // Prepare gex reference (Normal Ref)
        //
        if ( !cellranger_gex_index ) {

            // Make reference genome
            CELLRANGER_MKREF(
                ch_fasta,
                CELLRANGER_MKGTF.out.gtf,
                "gex_reference"
            )
            ch_versions = ch_versions.mix(CELLRANGER_MKREF.out.versions)
            ch_cellranger_gex_index = CELLRANGER_MKREF.out.reference.ifEmpty { empty_file }

        } else {
            ch_cellranger_gex_index = cellranger_gex_index
        }

        //
        // Prepare vdj reference (Special)
        //
        if ( !cellranger_vdj_index ) {

            if ( !params.skip_cellrangermulti_vdjref  ) { // if user uses cellranger multi but does not have VDJ data
                // Make reference genome
                CELLRANGER_MKVDJREF(
                    ch_fasta,
                    CELLRANGER_MKGTF.out.gtf,
                    [], // currently ignoring the 'seqs' option
                    "vdj_reference"
                )
                ch_versions = ch_versions.mix(CELLRANGER_MKVDJREF.out.versions)
                ch_cellranger_vdj_index = CELLRANGER_MKVDJREF.out.reference.ifEmpty { empty_file }
            } else {
                ch_cellranger_vdj_index = []
            }

        } else {
            ch_cellranger_vdj_index = cellranger_vdj_index
        }

        //
        // MODULE: cellranger multi
        //
        CELLRANGER_MULTI(
            ch_grouped_fastq.gex.map{ it[0] },
            ch_grouped_fastq.gex,
            ch_grouped_fastq.vdj,
            ch_grouped_fastq.ab,
            ch_grouped_fastq.beam,
            ch_grouped_fastq.cmo,
            ch_grouped_fastq.crispr,
            ch_cellranger_gex_index,
            ch_gex_frna_probeset,
            ch_gex_target_panel,
            ch_cellranger_vdj_index,
            ch_vdj_primer_index,
            ch_fb_reference,
            ch_beam_antigen_panel_csv,
            ch_beam_control_panel_csv,
            ch_gex_cmo_set,
            ch_cmo_barcode_csv,
            ch_gex_barcodes,
            ch_frna_sample_csv,
            params.skip_cellranger_renaming
        )
        ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.versions)

        //
        // Cellranger multi splits the results from each sample. So, a module execution will have: (1) a raw counts dir for all;
        // (2) a filtered counts dir PER sample.
        //
        // Thus, cellranger multi outputs data from all identified samples in a single channel, which will cause file collision.
        //
        // For the conversion, we should convert the resulting files of each sample, thus, now, we must parse the names
        // of the filtered 'per_sample_outs' of cellranger/multi and then, split the channels as such.
        //
        CELLRANGER_MULTI.out.outs
        .map{ meta, mtx_files ->
            def desired_files = []
            mtx_files.each{
                if ( it.toString().contains("per_sample_outs") ) {
                    def demuxed_sample_name = it.toString().split('per_sample_outs/')[1].split('/')[0]
                    def meta_clone    = meta.clone()
                    meta_clone.bkp_id = meta_clone.id
                    meta_clone.id     = demuxed_sample_name
                    desired_files.add( [ meta_clone, it ] )
                }
            }
            desired_files
        }
        .flatten()
        .buffer( size: 2 )   // gets: [ meta_clone, single_file ]
        .groupTuple( by: 0 ) // gets: [ meta_clone, all_files   ]
        .set { ch_parsed_per_sample_matrices }

        //
        // Split channels of raw and filtered to avoid file collision problems when loading the inputs in conversion modules.
        //
        ch_matrices_raw =
        CELLRANGER_MULTI.out.outs.map { meta, mtx_files -> // Use the unparsed matrix channel because raw results will be outside the 'per_sample_outs'
            def desired_files = []
            mtx_files.each{
                if (
                        it.toString().contains("raw_feature_bc_matrix") &&
                        // frna analysis also produces raw data per sample. For standardization, let's keep raw conversions always from the 'unmultiplexed' results.
                        !(it.toString().contains("sample_raw_feature_bc_matrix")) // TODO: Define if this is enough or if one still wants the raw conversion per samples when frna.
                    ) { desired_files.add( it ) }
            }
            [ meta, desired_files ]
        }

        ch_matrices_filtered =
        ch_parsed_per_sample_matrices.map { meta, mtx_files -> // Use the parsed matrix channel
            def desired_files = []
            mtx_files.each{
                if ( it.toString().contains("filtered_feature_bc_matrix") ) { desired_files.add( it ) }
            }
            [ meta, desired_files ]
        }

    emit:
        ch_versions
        cellrangermulti_out = CELLRANGER_MULTI.out.outs
        cellrangermulti_mtx = ch_matrices_raw.mix( ch_matrices_filtered )
}
