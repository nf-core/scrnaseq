/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { MTX_TO_H5AD   }             from '../../modules/local/mtx_to_h5ad.nf'
include { CONCAT_H5AD   }             from '../../modules/local/concat_h5ad.nf'
include { MTX_TO_SEURAT }             from '../../modules/local/mtx_to_seurat.nf'

workflow MTX_CONVERSION {

    take:
    mtx_matrices
    samplesheet

    main:
    //
    // Convert matrix do h5ad
    //
    MTX_TO_H5AD (
        mtx_matrices
    )

    //
    // Concat sample-specific h5ad in one
    //
    CONCAT_H5AD (
        MTX_TO_H5AD.out.h5ad.collect(), // gather all sample-specific files
        samplesheet
    )

    //
    // Convert matrix do seurat
    //
    MTX_TO_SEURAT (
        mtx_matrices
    )

    emit:
    MTX_TO_SEURAT.out.versions

}
