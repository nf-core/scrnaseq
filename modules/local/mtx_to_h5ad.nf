process MTX_TO_H5AD {
    tag "$prefix"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gcfntnu/scanpy:1.7.0' :
        'gcfntnu/scanpy:1.7.0' }"

    input:
    tuple val(cellranger_prefix), path(cellranger_outdir)

    output:
    path "matrix.h5ad", emit: h5ad

    script:
    def prefix = cellranger_prefix.tokenize('-')[1]
    """
    mkdir -p ${cellranger_prefix}
    mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o ${cellranger_prefix}/matrix.h5ad
    """
}
