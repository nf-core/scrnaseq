/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { MTX_TO_H5AD   }             from '../../modules/local/mtx_to_h5ad.nf'
include { CONCAT_H5AD   }             from '../../modules/local/concat_h5ad.nf'
include { MTX_TO_SEURAT }             from '../../modules/local/mtx_to_seurat.nf'

workflow MTX_CONVERSION {

    take:
    mtx_matrices
    samplesheet
    txp2gene
    star_index

    main:
        ch_versions = Channel.empty()

        // Cellranger module output contains too many files which cause path collisions, we filter to the ones we need.
        if ( params.aligner == "cellranger" ||  params.aligner == "cellrangerarc" ) {
            mtx_matrices = mtx_matrices.map { meta, mtx_files ->
                    [ meta, mtx_files.findAll { it.toString().contains("filtered_feature_bc_matrix") } ]
                }
                .filter { meta, mtx_files -> mtx_files } // Remove any that are missing the relevant files
        }

        //
        // Convert matrix to h5ad
        //
        MTX_TO_H5AD (
            mtx_matrices,
            txp2gene,
            star_index
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

        //TODO CONCAT h5ad and MTX to h5ad should also have versions.yaml output
        ch_version = ch_versions.mix(MTX_TO_H5AD.out.versions, MTX_TO_SEURAT.out.versions)

    emit:
    ch_versions
    counts = MTX_TO_H5AD.out.counts

}
