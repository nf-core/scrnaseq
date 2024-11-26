/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { CONCAT_H5AD           } from '../../modules/local/concat_h5ad.nf'
include { ANNDATAR_CONVERT      } from '../../modules/local/anndatar_convert'

workflow H5AD_CONVERSION {

    take:
    ch_h5ads
    samplesheet

    main:
        ch_versions = Channel.empty()

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
        ch_versions = ch_versions.mix(CONCAT_H5AD.out.versions.first())

        //
        // MODULE: Convert to Rds with AnndataR package
        //
        ANNDATAR_CONVERT (
            ch_h5ads.mix( ch_h5ad_concat )
        )
        ch_versions = ch_versions.mix(ANNDATAR_CONVERT.out.versions.first())

    emit:
    ch_versions
    h5ads = ch_h5ads

}
