process MTX_TO_H5AD {
    tag "$prefix"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gcfntnu/scanpy:1.7.0' :
        'gcfntnu/scanpy:1.7.0' }"

    input:
    path cellranger_outdir

    output:
    path "matrix.h5ad", emit: h5ad

    script:
    def prefix = cellranger_outdir[0].getName().toString()[-3]
    """
    mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o matrix.h5ad
    """
}
