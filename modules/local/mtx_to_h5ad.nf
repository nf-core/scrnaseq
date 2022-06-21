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
    if (params.aligner == 'cellranger')
    """
    # convert file types
    cellranger_mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o matrix.h5ad
    """

    else if (params.aligner == 'kallisto')
    """
    # convert file types
    mtx_to_h5ad.py \\
        -m *_kallistobustools_count/counts_unfiltered/*.mtx \\
        -b *_kallistobustools_count/counts_unfiltered/*.barcodes.txt \\
        -f *_kallistobustools_count/counts_unfiltered/*.genes.txt \\
        -o cells_x_genes.h5ad
    """

    else if (params.aligner == 'alevin')
    """
    # convert file types
    mtx_to_h5ad.py \\
        -m *_alevin_results/alevin/quants_mat.mtx.gz \\
        -b *_alevin_results/alevin/quants_mat_rows.txt \\
        -f *_alevin_results/alevin/quants_mat_cols.txt \\
        -o quants_mat.h5ad
    """

    stub:
    """
    touch matrix.h5ad
    """
}
