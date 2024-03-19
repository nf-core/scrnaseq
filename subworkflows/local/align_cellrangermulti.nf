//
// Include modules
//
include { CELLRANGER_MKGTF    } from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF    } from "../../modules/nf-core/cellranger/mkref/main.nf"
include { CELLRANGER_MKVDJREF } from "../../modules/nf-core/cellranger/mkvdjref/main.nf"
include { CELLRANGER_MULTI    } from "../../modules/nf-core/cellranger/multi/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI_ALIGN {
    take:
        ch_fasta
        ch_gtf
        ch_fastq
        cellranger_gex_index
        cellranger_vdj_index
        empty_file

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
        ch_beam_antigen_panel_csv = params.beam_antigen_csv   ? file(params.beam_antigen_csv)   : empty_file
        ch_beam_control_panel_csv = params.beam_control_csv   ? file(params.beam_control_csv)   : empty_file
        ch_frna_sample_csv        = params.frna_sample_csv    ? file(params.frna_sample_csv)    : empty_file
        ch_cmo_barcode_csv        = params.cmo_barcode_csv    ? file(params.cmo_barcode_csv)    : empty_file

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

            // Make reference genome
            CELLRANGER_MKVDJREF(
                ch_fasta,
                CELLRANGER_MKGTF.out.gtf,
                "vdj_reference"
            )
            ch_versions = ch_versions.mix(CELLRANGER_MKVDJREF.out.versions)
            ch_cellranger_vdj_index = CELLRANGER_MKVDJREF.out.reference.ifEmpty { empty_file }

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

    emit:
        ch_versions
        cellrangermulti_out = CELLRANGER_MULTI.out.outs
}
