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
    path "matrix.h5ad.gz", emit: h5ad

    script:
    if (params.aligner == 'cellranger')
    """
    # convert file types
    cellranger_mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o matrix.h5ad

    gzip -c matrix.h5ad > matrix.h5ad.gz
    """

    else if (params.aligner == 'kallisto')
    """
    # convert file types
    mtx_to_h5ad.py \\
        -m ${meta.id}_kallistobustools_count/counts_unfiltered/*.mtx \\
        -b ${meta.id}_kallistobustools_count/counts_unfiltered/*.barcodes.txt \\
        -f ${meta.id}_kallistobustools_count/counts_unfiltered/*.genes.txt \\
        -o matrix.h5ad

    gzip -c matrix.h5ad > matrix.h5ad.gz
    """

    stub:
    """
    touch matrix.h5ad.gz
    """
}
