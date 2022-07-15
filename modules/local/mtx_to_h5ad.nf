process MTX_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)

    output:
    path "*.h5ad", emit: h5ad

    script:
    // def file paths for aligners (except cellranger)
    if (params.aligner == 'kallisto') {
        mtx_matrix   = "*count/counts_unfiltered/*.mtx"
        barcodes_tsv = "*count/counts_unfiltered/*.barcodes.txt"
        features_tsv = "*count/counts_unfiltered/*.genes.txt"
    } else if (params.aligner == 'alevin') {
        mtx_matrix   = "*_alevin_results/alevin/quants_mat.mtx.gz"
        barcodes_tsv = "*_alevin_results/alevin/quants_mat_rows.txt"
        features_tsv = "*_alevin_results/alevin/quants_mat_cols.txt"
    } else if (params.aligner == 'star') {
        mtx_matrix   = "*.Solo.out/Gene/filtered/matrix.mtx"
        barcodes_tsv = "*.Solo.out/Gene/filtered/barcodes.tsv"
        features_tsv = "*.Solo.out/Gene/filtered/features.tsv"
    }

    //
    // run script
    //
    if (params.aligner == 'cellranger')
    """
    # convert file types
    cellranger_mtx_to_h5ad.py \\
        --mtx filtered_feature_bc_matrix.h5 \\
        --sample ${meta.id} \\
        --out ${meta.id}_matrix.h5ad
    """

    else
    """
    # convert file types
    mtx_to_h5ad.py \\
        --aligner ${params.aligner} \\
        --sample ${meta.id} \\
        --mtx $mtx_matrix \\
        --barcode $barcodes_tsv \\
        --feature $features_tsv \\
        --out ${meta.id}_matrix.h5ad
    """

    stub:
    """
    touch ${meta.id}_matrix.h5ad
    """
}
