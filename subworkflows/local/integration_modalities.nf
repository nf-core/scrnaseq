/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { INTEGRATION } from '../../modules/local/integration'
include { MOFA_INTEGRATION } from '../../modules/local/MOFA_integration'

workflow INTEGRATION_MODALITIES {
    take:
    h5mu
    h5ad
    n_neighbors_harmony
    min_dist_harmony
    integration_var

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Integrate GEX and ADT modalities indipendently with Harmony
    //
    INTEGRATION(
        h5mu,
        params.n_neighbors_harmony,
        params.min_dist_harmony,
        params.integration_var,
    )
    ch_versions = ch_versions.mix(INTEGRATION.out.versions.first())
    integration_out = INTEGRATION.out.h5mu
    //
    // MODULE: Compute WNN graph for each modality
    //
    if (!params.skip_integration) {
        MOFA_INTEGRATION(
            INTEGRATION.out.h5mu,
            h5ad,
        )
        ch_versions = ch_versions.mix(MOFA_INTEGRATION.out.versions.first())
        mofa_out = MOFA_INTEGRATION.out.h5mu
        h5mu_out = integration_out.mix(mofa_out)
    }
    else {
        h5mu_out = integration_out
    }

    emit:
    ch_versions
    h5mu_out
}
