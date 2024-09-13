/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { ANNDATAR_CONVERT      } from '../../modules/local/anndatar_convert'
include { MTX_TO_H5AD           } from '../../modules/local/mtx_to_h5ad.nf'
include { CONCAT_H5AD           } from '../../modules/local/concat_h5ad.nf'
include { MTX_TO_SEURAT         } from '../../modules/local/mtx_to_seurat.nf'
include { EMPTY_DROPLET_REMOVAL } from './emptydrops_removal'

workflow MTX_CONVERSION {

    take:
    mtx_matrices
    samplesheet

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Standardize h5ad and convert with AnndataR package
        //
        ANNDATAR_CONVERT (
            mtx_matrices
        )

        //
        // Concat sample-specific h5ad in one
        //
        ch_concat_h5ad_input = mtx_matrices.groupTuple() // gather all sample-specific files / per type
        if (params.aligner == 'kallisto' && params.kb_workflow != 'standard') {
            // when having spliced / unspliced matrices, the collected tuple has two levels ( [[mtx_1, mtx_2]] )
            // which nextflow break because it is not a valid 'path' thus, we have to remove one level
            // making it as [ mtx_1, mtx_2 ]
            ch_concat_h5ad_input = ch_concat_h5ad_input.map{ type, matrices -> [ type, matrices.flatten().toList() ] }
        }

        CONCAT_H5AD (
            ch_concat_h5ad_input,
            samplesheet
        )

        //TODO CONCAT h5ad and MTX to h5ad should also have versions.yaml output
        // ch_versions = ch_versions.mix(MTX_TO_H5AD.out.versions, MTX_TO_SEURAT.out.versions)

    emit:
    ch_versions
    // counts = MTX_TO_H5AD.out.counts  was this ever used?

}
