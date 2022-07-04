process MTX_TO_SEURAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "seurat-scripts" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://satijalab/seurat:4.1.0' :
        'satijalab/seurat:4.1.0' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)

    output:
    path "*.seurat", emit: h5ad

    script:
    def aligner = params.aligner
    if (params.aligner == "cellranger") {
        matrix   = "filtered_feature_bc_matrix/matrix.mtx.gz"
        barcodes = "filtered_feature_bc_matrix/barcodes.tsv.gz"
        features = "filtered_feature_bc_matrix/features.tsv.gz"
    } else if (params.aligner == "kallisto") {
        matrix   = "*_kallistobustools_count/counts_unfiltered/*.mtx"
        barcodes = "*_kallistobustools_count/counts_unfiltered/*.barcodes.txt"
        features = "*_kallistobustools_count/counts_unfiltered/*.genes.txt"
    } else if (params.aligner == "alevin") {
        matrix   = "*_alevin_results/alevin/quants_mat.mtx.gz"
        barcodes = "*_alevin_results/alevin/quants_mat_rows.txt"
        features = "*_alevin_results/alevin/quants_mat_cols.txt"
    } else if (params.aligner == 'star') {
        matrix   = "*.Solo.out/Gene/filtered/matrix.mtx"
        barcodes = "*.Solo.out/Gene/filtered/barcodes.tsv"
        features = "*.Solo.out/Gene/filtered/features.tsv"
    }

    """
    mtx_to_seurat.R \\
        $matrix \\
        $barcodes \\
        $features \\
        ${meta.id}_matrix.seurat \\
        ${aligner}
    """

    stub:
    """
    touch ${meta.id}_matrix.seurat
    """
}
