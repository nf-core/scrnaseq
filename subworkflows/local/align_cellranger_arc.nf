/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGER_ARC_MKGTF} from "../../modules/local/cellranger_arc/mkgtf/main.nf"
include {CELLRANGER_ARC_MKREF} from "../../modules/local/cellranger_arc/mkref/main.nf"
include {GENERATE_LIB_CSV} from "../../modules/local/generate_cellranger_lib_csv.nf"
include {CELLRANGER_ARC_COUNT} from "../../modules/local/cellranger_arc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ARC_ALIGN {
    take:
        fasta
        gtf
        motifs
        reference_config
        cellranger_index
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf && motifs && reference_config):
            "Must provide either a cellranger-atac index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf') + motif file (--motifs) + reference_config ('--reference_config')."

        if (!cellranger_arc_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_ARC_MKGTF( gtf )
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKGTF.out.versions)

            // Make reference genome
            CELLRANGER_ARC_MKREF( fasta, CELLRANGER_ARC_MKGTF.out.gtf, motifs, reference_config, "cellranger_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKREF.out.versions)
            cellranger_arc_index = CELLRANGER_ARC_MKREF.out.reference
        }

        //TOFLO do I need to copy the files just for the meta data?
        GENERATE_LIB_CSV( ch_fastq )

        // Obtain read counts
        CELLRANGER_ARC_COUNT (
            ch_fastq,
            GENERATE_LIB_CSV.out.csv,
            cellranger_arc_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_ARC_COUNT.out.versions)
        
    emit:
        ch_versions
        cellranger_arc_out  = CELLRANGER_ARC_COUNT.out.outs
}