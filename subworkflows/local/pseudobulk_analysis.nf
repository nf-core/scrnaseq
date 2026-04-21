include { PSEUDOBULK } from '../../modules/local/pseudobulk/main.nf'
include { DIFFERENTIAL_ANALYSIS } from '../../subworkflows/nfdata-omics/deseq2_analysis/main.nf'

workflow PSEUDOBULK_ANALYSIS {
    take:
    h5mu
    resolution
    group_column
    comparisons
    formula
    fdr

    main:

    ch_versions = channel.empty()

    h5mu
        .combine(resolution)
        .combine(group_column)
        .combine(comparisons)
        .map { meta, file, res, column, comp -> tuple(meta, file, column, res, comp) }
        .set { pseudobulk_inputs }

    PSEUDOBULK(
        pseudobulk_inputs
    )

    ch_versions = ch_versions.mix(PSEUDOBULK.out.versions)

    PSEUDOBULK.out.pseudobulk_deseq2
        .flatMap { meta, export_dir ->
            export_dir
                .listFiles()
                .collect { cluster_dir ->
                    def cluster = cluster_dir.name
                    def coldata = cluster_dir.resolve("coldata_cl_${cluster}.tsv")
                    def counts = cluster_dir.resolve("counts_cl_${cluster}.tsv")
                    [[id: cluster, resolution: meta], counts, coldata]
                }
        }
        .multiMap { cluster, counts, coldata ->
            ch_counts: [cluster, counts]
            ch_metadata: [cluster, coldata]
        }
        .set { deseq2_inputs }

    DIFFERENTIAL_ANALYSIS(
        deseq2_inputs.ch_counts,
        deseq2_inputs.ch_metadata,
        formula,
        comparisons,
        fdr,
    )

    ch_versions = ch_versions.mix(DIFFERENTIAL_ANALYSIS.out.versions)

    emit:
    pseudobulk_plots = PSEUDOBULK.out.pseudobulk_plots
    pseudobulk_tables = PSEUDOBULK.out.pseudobulk_deseq2
    rds = DIFFERENTIAL_ANALYSIS.out.rds
    dge = DIFFERENTIAL_ANALYSIS.out.dge
    summary = DIFFERENTIAL_ANALYSIS.out.summary
    hist_pval = DIFFERENTIAL_ANALYSIS.out.hist_pval
    volcano = DIFFERENTIAL_ANALYSIS.out.volcano
    versions = ch_versions
}
