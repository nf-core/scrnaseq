//
// Include modules
//
include { CELLRANGER_MKGTF                  } from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF                  } from "../../modules/nf-core/cellranger/mkref/main.nf"
include { CELLRANGER_MKVDJREF               } from "../../modules/nf-core/cellranger/mkvdjref/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI_REF {
    take:
        ch_fasta
        ch_gtf
        cellranger_gex_index
        cellranger_vdj_index

    main:
        ch_versions    = Channel.empty()

        //
        // TODO: Include checkers for cellranger multi parameter combinations. For example, when VDJ data is given, require VDJ ref. If FFPE, require frna probe sets, etc.
        //

        //
        // Prepare GTF
        //
        if ( !cellranger_gex_index || (!cellranger_vdj_index && !params.skip_cellrangermulti_vdjref) ) {

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
            ch_cellranger_gex_index = CELLRANGER_MKREF.out.reference.ifEmpty { [] }

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
                ch_cellranger_vdj_index = CELLRANGER_MKVDJREF.out.reference.ifEmpty { [] }
            } else {
                ch_cellranger_vdj_index = []
            }

        } else {
            ch_cellranger_vdj_index = cellranger_vdj_index
        }

    emit:
        ch_versions
        ch_cellranger_gex_index
        ch_cellranger_vdj_index
}
