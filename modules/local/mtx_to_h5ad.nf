process MTX_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gcfntnu/scanpy:1.7.0' :
        'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0' }"

    input:
    tuple val(meta), path(cellranger_outdir)

    output:
    path "${meta.cellranger_prefix}", emit: h5ad

    script:
    """
    mkdir -p ${meta.cellranger_prefix}
    mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o ${meta.cellranger_prefix}/matrix.h5ad
    """
}
