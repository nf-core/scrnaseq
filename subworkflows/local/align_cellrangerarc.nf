/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGERARC_MKGTF} from "../../modules/nf-core/cellrangerarc/mkgtf/main.nf"
include {CELLRANGERARC_MKREF} from "../../modules/nf-core/cellrangerarc/mkref/main.nf"
include {CELLRANGERARC_COUNT} from "../../modules/nf-core/cellrangerarc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGERARC_ALIGN {
    take:
        fasta
        gtf
        motifs
        cellranger_index
        ch_fastq
        cellrangerarc_config

    main:
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf')."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGERARC_MKGTF( gtf )
            filtered_gtf = CELLRANGERARC_MKGTF.out.gtf
            ch_versions = ch_versions.mix(CELLRANGERARC_MKGTF.out.versions)

            // Make reference genome
            assert (( !params.cellrangerarc_reference && !cellrangerarc_config ) ||
                    ( params.cellrangerarc_reference && cellrangerarc_config ) ) :
                "If you provide a config file you also have to specific the reference name and vice versa."

            cellrangerarc_reference = 'cellrangerarc_reference'
            if ( params.cellrangerarc_reference ){
                cellrangerarc_reference = params.cellrangerarc_reference
            }

            CELLRANGERARC_MKREF( fasta, filtered_gtf, motifs, cellrangerarc_config, cellrangerarc_reference )
            ch_versions = ch_versions.mix(CELLRANGERARC_MKREF.out.versions)
            cellranger_index = CELLRANGERARC_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGERARC_COUNT (
            ch_fastq,
            cellranger_index
        )
        ch_versions = ch_versions.mix(CELLRANGERARC_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_arc_out  = CELLRANGERARC_COUNT.out.outs
}
