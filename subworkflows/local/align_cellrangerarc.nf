/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGERARC_MKGTF} from "../../modules/local/cellrangerarc/mkgtf/main.nf"
include {CELLRANGERARC_MKREF} from "../../modules/local/cellrangerarc/mkref/main.nf"
include {CELLRANGERARC_GENERATECONFIG} from "../../modules/local/generate_cellranger_mkref_config.nf"
include {CELLRANGERARC_COUNT} from "../../modules/local/cellrangerarc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGERARC_ALIGN {
    take:
        fasta
        gtf
        motifs
        cellrangerarc_index
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert cellrangerarc_index || (fasta && gtf && motifs):
            "Must provide either a cellranger-atac index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf') + motif file (--motifs)."

        if (!cellrangerarc_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGERARC_MKGTF( gtf )
            filtered_gtf = CELLRANGERARC_MKGTF.out.gtf
            ch_versions = ch_versions.mix(CELLRANGERARC_MKGTF.out.versions)

            // Generate the config for mkref
            CELLRANGERARC_GENERATECONFIG(fasta.name, filtered_gtf.name, motifs.name)
            ch_versions = ch_versions.mix(CELLRANGERARC_GENERATECONFIG.out.versions)

            // Make reference genome
            CELLRANGERARC_MKREF( fasta, filtered_gtf, motifs, CELLRANGERARC_GENERATECONFIG.out.config, "cellrangerarc_reference" )
            ch_versions = ch_versions.mix(CELLRANGERARC_MKREF.out.versions)
            cellrangerarc_index = CELLRANGERARC_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGERARC_COUNT (
            ch_fastq,
            cellrangerarc_index
        )
        ch_versions = ch_versions.mix(CELLRANGERARC_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_arc_out  = CELLRANGERARC_COUNT.out.outs
}