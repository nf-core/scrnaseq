/*
 * Alignment with Cellranger
 */

include {UNIVERSC_MKGTF} from "../../modules/nf-core/universc/mkgtf/main.nf"
include {UNIVERSC_MKREF} from "../../modules/nf-core/universc/mkref/main.nf"
include {UNIVERSC_LAUNCH} from "../../modules/nf-core/universc/launch/main.nf"
include {MTX_TO_H5AD     } from "../../modules/local/mtx_to_h5ad.nf"

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
            UNIVERSC_MKGTF( gtf )
            ch_versions = ch_versions.mix(UNIVERSC_MKGTF.out.versions)

            // Make reference genome
            UNIVERSC_MKREF( fasta, UNIVERSC_MKGTF.out.gtf, "cellranger_reference" )
            ch_versions = ch_versions.mix(UNIVERSC_MKREF.out.versions)
            universc_index = UNIVERSC_MKREF.out.reference
        }

        // Obtain read counts
        UNIVERSC_LAUNCH (
            // TODO add technology and chemistry input parameters and set defaults
            ch_fastq.map{ meta, reads -> [meta + ["id": meta.id, "samples:", meta.id, "technology": universc_technology, "single_end:", false, "strandedness:", 'forward'], reads] },
            universc_index
        )
        ch_versions = ch_versions.mix(UNIVERSC_LAUNCH.out.versions)

    emit:
        ch_versions
        universc_out  = UNIVERSC_LAUNCH.out.outs
}
