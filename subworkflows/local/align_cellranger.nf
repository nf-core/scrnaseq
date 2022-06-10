/*
 * Alignment with Cellranger
 */

include {CELLRANGER_MKGTF} from "../../modules/nf-core/modules/cellranger/mkgtf/main.nf"
include {CELLRANGER_MKREF} from "../../modules/nf-core/modules/cellranger/mkref/main.nf"
include {CELLRANGER_COUNT} from "../../modules/nf-core/modules/cellranger/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ALIGN {
    take:
        fasta
        gtf
        cellranger_index
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--genome_fasta') and a gtf file ('--gtf')."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF( gtf )
            ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

            CELLRANGER_MKGTF.out.gtf.view()
            // Make reference genome
            CELLRANGER_MKREF( fasta, CELLRANGER_MKGTF.out.gtf, "cellranger_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_MKREF.out.versions)
            cellranger_index = CELLRANGER_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGER_COUNT (
            // TODO what is `gem` and why is it needed?
            ch_fastq.map{ meta, reads -> [meta + ["gem": meta.id, "samples": [meta.id]], reads] },
            cellranger_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_out  = CELLRANGER_COUNT.out.outs
}
