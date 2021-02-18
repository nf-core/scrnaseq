/*
 * Alignment with Cellranger
 */

params.cellranger_mkgtf_options    = [:]
params.cellranger_mkref_options    = [:]
params.cellranger_count_options    = [:]

include {CELLRANGER_MKGTF} from "../../modules/local/software/cellranger/mkgtf/main.nf" addParams(options: params.cellranger_mkgtf_options)
include {CELLRANGER_MKREF} from "../../modules/local/software/cellranger/mkref/main.nf" addParams(options: params.cellranger_mkref_options)
include {CELLRANGER_COUNT} from "../../modules/local/software/cellranger/count/main.nf" addParams(options: params.cellranger_count_options)

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ALIGN {
    take:
        fasta
        gtf
        reads

    main:
    
        // Filter GTF based on gene biotypes passed in params.modules
        CELLRANGER_MKGTF( gtf )

        // Make reference genome
        CELLRANGER_MKREF( CELLRANGER_MKGTF.out, fasta )

        // Obtain read counts
        CELLRANGER_COUNT( reads, CELLRANGER_MKREF.out.collect() )

    emit:
        read_counts     = CELLRANGER_COUNT.out.read_counts
        cellranger_out  = CELLRANGER_COUNT.out.cellranger_out
}