/*
 * Alignment with Cellranger open-source implementation called by UniverSC
 */

include {CELLRANGER_MKGTF} from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include {CELLRANGER_MKREF} from "../../modules/nf-core/cellranger/mkref/main.nf"
include {UNIVERSC} from "../../modules/nf-core/universc/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow UNIVERSC_ALIGN {
    take:
        fasta
        gtf
        universc_index
        universc_technology
        ch_fastq

    main:
        ch_versions = Channel.empty()

        assert universc_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--fasta') and a gtf file ('--gtf')."

        if (!universc_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF( gtf )
            ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

            // Make reference genome
            CELLRANGER_MKREF( fasta, CELLRANGER_MKGTF.out.gtf, "cellranger_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_MKREF.out.versions)
            universc_index = CELLRANGER_MKREF.out.reference
        }

        // Obtain read counts
        UNIVERSC (
            ch_fastq.map{
                meta, reads -> [
                    // defaults
                    ["samples": [meta.id], "technology": universc_technology, "chemistry": "auto", "single_end": false, "strandedness": "forward"] + meta, // + meta overrides defaults with information already in meta
                    reads
                ]
            },
            universc_index
        )
        ch_versions = ch_versions.mix(UNIVERSC.out.versions)

    emit:
        ch_versions
        universc_out  = UNIVERSC.out.outs
        star_index = universc_index
}
