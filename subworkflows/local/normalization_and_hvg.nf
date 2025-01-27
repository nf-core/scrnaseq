/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { NORMALIZATION } from '../../modules/local/normalization'
include { HIGHLY_VARIABLE_GENES } from '../../modules/local/highly_variable_genes'

workflow H5AD_CONVERSION {

    take:
    ch_h5ads

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Normalize count matrices contained in the concatenated h5ad file
        //
        NORMALIZATION (
            ch_h5ads
        )
        ch_versions = ch_versions.mix(NORMALIZATION.out.versions.first())

        //
        // MODULE: Highly variable genes detection, added with gene annotation
        //
        HIGHLY_VARIABLE_GENES (
            NORMALIZATION.out.h5ad
        )
        ch_versions = ch_versions.mix(HIGHLY_VARIABLE_GENES.out.versions.first())

    emit:
    ch_versions
    h5ads = HIGHLY_VARIABLE_GENES.out.h5ad

}
