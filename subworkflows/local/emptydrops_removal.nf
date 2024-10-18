include { CELLBENDER_REMOVEBACKGROUND } from '../../modules/nf-core/cellbender/removebackground'
include { ADATA_BARCODES              } from '../../modules/local/adata_barcodes'

workflow EMPTY_DROPLET_REMOVAL {
    take:
    ch_unfiltered

    main:
    ch_versions = Channel.empty()

    CELLBENDER_REMOVEBACKGROUND(ch_unfiltered)
    ch_versions = ch_versions.mix(CELLBENDER_REMOVEBACKGROUND.out.versions)

    ch_combined =
    ch_unfiltered
        .join(CELLBENDER_REMOVEBACKGROUND.out.barcodes)
        .map { meta, h5ad, csv ->
            def meta_clone = meta.clone()
            meta_clone.input_type = meta['input_type'].toString().replaceAll('raw', 'emptydrops_filter')

            [ meta_clone, h5ad, csv ]
        }

    ADATA_BARCODES(ch_combined)
    ch_versions = ch_versions.mix(ADATA_BARCODES.out.versions)

    ch_h5ad = ADATA_BARCODES.out.h5ad

    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}
