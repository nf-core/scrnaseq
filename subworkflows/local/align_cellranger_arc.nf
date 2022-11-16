/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGER_ARC_MKGTF} from "../../modules/local/cellranger_arc/mkgtf/main.nf"
include {CELLRANGER_ARC_MKREF} from "../../modules/local/cellranger_arc/mkref/main.nf"
include {GENERATE_LIB_CSV} from "../../modules/local/generate_cellranger_lib_csv.nf"
include {GENERATE_CONFIG} from "../../modules/local/generate_cellranger_mkref_config.nf"
include {CELLRANGER_ARC_COUNT} from "../../modules/local/cellranger_arc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ARC_ALIGN {
    take:
        fasta
        gtf
        motifs
        cellranger_index
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf && motifs):
            "Must provide either a cellranger-atac index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf') + motif file (--motifs)."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_ARC_MKGTF( gtf )
            filtered_gtf = CELLRANGER_ARC_MKGTF.out.gtf
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKGTF.out.versions)

            // Generate the config for mkref
            GENERATE_CONFIG(fasta.name, filtered_gtf.name, motifs.name)
            ch_versions.mix(GENERATE_CONFIG.out.versions)

            // Make reference genome
            CELLRANGER_ARC_MKREF( fasta, filtered_gtf, motifs, GENERATE_CONFIG.out.config, "cellranger_arc_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKREF.out.versions)
            cellranger_index = CELLRANGER_ARC_MKREF.out.reference
        }

        //TOFLO do I need to copy the files just for the meta data?
        //GENERATE_LIB_CSV( ch_fastq )
        //ch_versions.mix(GENERATE_LIB_CSV.out.versions)

        // Obtain read counts
        CELLRANGER_ARC_COUNT (
            ch_fastq,
            cellranger_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_ARC_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_arc_out  = CELLRANGER_ARC_COUNT.out.outs
}