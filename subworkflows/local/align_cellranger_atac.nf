/*
 * Alignment with Cellranger ATAC
 */
include {CELLRANGER_ATAC_MKREF} from "../../modules/local/cellranger_atac/mkref/main.nf"
include {CELLRANGER_ATAC_COUNT} from "../../modules/local/cellranger_atac/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ATAC_ALIGN {
    take:
        reference_config
        cellranger_atac_index
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert cellranger_atac_index || reference_config:
            "Must provide either a cellranger-atac index or reference_config ('--reference_config')."

        if (!cellranger_atac_index) {
            // Make reference genome
            CELLRANGER_ATAC_MKREF( reference_config, "cellranger_atac_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_ATAC_MKREF.out.versions)
            cellranger_atac_index = CELLRANGER_ATAC_MKREF.out.reference
        }

        /*
        // Obtain read counts
        CELLRANGER_ATAC_COUNT (
            // TODO what is `gem` and why is it needed?
            ch_fastq.map{ meta, reads -> [meta + ["gem": meta.id, "samples": [meta.id]], reads] },
            cellranger_atac_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_ATAC_COUNT.out.versions)
        */


    emit:
        ch_versions
        cellranger_atac_index
        //cellranger_atac_out  = CELLRANGER_ATAC_COUNT.out.outs
}