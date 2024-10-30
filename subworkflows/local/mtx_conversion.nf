/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { MTX_TO_H5AD           } from '../../modules/local/mtx_to_h5ad'
include { CONCAT_H5AD           } from '../../modules/local/concat_h5ad.nf'
include { ANNDATAR_CONVERT      } from '../../modules/local/anndatar_convert'
include { EMPTY_DROPLET_REMOVAL } from '../../subworkflows/local/emptydrops_removal'

workflow MTX_CONVERSION {

    take:
    mtx_matrices
    txp2gene
    star_index
    samplesheet

    main:
        ch_versions = Channel.empty()
        ch_h5ads    = Channel.empty()

        //
        // MODULE: Convert matrices to h5ad
        //
        MTX_TO_H5AD (
            mtx_matrices,
            txp2gene,
            star_index
        )
        ch_versions = ch_versions.mix(MTX_TO_H5AD.out.versions.first())

        // fix channel size when kallisto non-standard workflow
        if (params.aligner == 'kallisto' && !(params.kb_workflow == 'standard')) {
            ch_h5ads =
            MTX_TO_H5AD.out.h5ad
                .transpose()
                .map { meta, h5ad ->
                    def meta_clone = meta.clone()
                    def spc_prefix = h5ad.toString().contains('unspliced') ? 'un' : ''

                    meta_clone["input_type"] = "${meta.input_type}_${spc_prefix}spliced"

                    [ meta_clone, h5ad ]
                }
        } else {
            ch_h5ads = MTX_TO_H5AD.out.h5ad
        }

        //
        // SUBWORKFLOW: Run cellbender emptydrops filter
        //
        if ( !params.skip_emptydrops && !(params.aligner in ['cellrangerarc']) ) {

            // emptydrops should only run on the raw matrices thus, filter-out the filtered result of the aligners that can produce it
            EMPTY_DROPLET_REMOVAL (
                ch_h5ads.filter { meta, mtx_files -> meta.input_type.contains('raw') }
            )
            ch_h5ads = ch_h5ads.mix( EMPTY_DROPLET_REMOVAL.out.h5ad )

        }

        //
        // Concat sample-specific h5ad in one
        //
        ch_concat_h5ad_input = ch_h5ads.groupTuple() // gather all sample-specific files / per type
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
        ch_h5ad_concat = CONCAT_H5AD.out.h5ad.map{ meta, file ->
            def meta_clone = meta.clone()
            meta_clone.id = 'combined' // maintain output prefix
            [ meta_clone, file ]
        }

        //
        // MODULE: Convert to Rds with AnndataR package
        //
        ANNDATAR_CONVERT (
            ch_h5ads.mix( ch_h5ad_concat )
        )

        //TODO CONCAT h5ad and MTX to h5ad should also have versions.yaml output
        // ch_versions = ch_versions.mix(MTX_TO_H5AD.out.versions, MTX_TO_SEURAT.out.versions)

    emit:
    ch_versions
    // counts = MTX_TO_H5AD.out.counts  was this ever used?

}
