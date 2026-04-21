include { DESEQ2_FIT } from "../../../modules/nfdata-omics/deseq2/fit/main.nf"
include { DESEQ2_COMPARE } from "../../../modules/nfdata-omics/deseq2/compare/main.nf"

workflow DIFFERENTIAL_ANALYSIS {
    take:
    ch_counts // channel: [val(meta), path(counts)]
    ch_metadata // channel: [val(meta), path(metadata)]
    ch_formula // channel: val(formula)
    ch_comparison // channel: val(comparison)
    ch_fdr_threshold // channel: val(fdr_threshold)

    main:

    ch_versions = channel.empty()

    DESEQ2_FIT(
        ch_counts,
        ch_metadata,
        ch_formula,
    )

    ch_versions = ch_versions.mix(DESEQ2_FIT.out.versions)

    ch_rds_comparison = DESEQ2_FIT.out.rds
        .combine(ch_comparison)
        .map { meta, rds, comparison ->
            tuple(meta, rds, comparison)
        }

    DESEQ2_COMPARE(
        ch_rds_comparison,
        ch_fdr_threshold,
    )

    ch_versions = ch_versions.mix(DESEQ2_COMPARE.out.versions)

    emit:
    rds = DESEQ2_FIT.out.rds // channel: [val(meta), path(deseq2_obj.rds)]
    dge = DESEQ2_COMPARE.out.dge // channel: [val(meta), path(deseq2_toptable.*.txt)]
    summary = DESEQ2_COMPARE.out.summary // channel: [val(meta), path(deseq2_summary.*.txt)]
    hist_pval = DESEQ2_COMPARE.out.hist_pval // channel: [val(meta), path(pvalue_hist.*.pdf)]
    volcano = DESEQ2_COMPARE.out.volcano // channel: [val(meta), path(volcano.*.pdf)]
    versions = ch_versions // channel: versions.yml
}
