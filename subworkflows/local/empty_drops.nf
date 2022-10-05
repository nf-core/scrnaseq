/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { EMPTYDROPS_FILTER   }             from '../../modules/local/empty_drops_filter.nf'

workflow emptydrops {

    take:
    mtx_matrices

    main:
    //
    // Filter matrix by emptydrops
    //
    EMPTYDROPS_FILTER (
        mtx_matrices
    )