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
    path "${meta.cellranger_prefix}/outs/filtered_feature_bc_matrix/matrix.h5ad.gz", emit: h5ad

    script:
    """
    # create dir to mirror cellranger output organisation to have results published in the same place
    mkdir -p ${meta.cellranger_prefix}/outs/filtered_feature_bc_matrix ;

    # convert file types
    mtx_to_h5ad.py \\
        -m filtered_feature_bc_matrix \\
        -o matrix.h5ad
    
    gzip -c matrix.h5ad > ${meta.cellranger_prefix}/outs/filtered_feature_bc_matrix/matrix.h5ad.gz
    """

    stub:
    """
    # create dir to mirror cellranger output organisation to have results published in the same place
    mkdir -p ${meta.cellranger_prefix}/outs/filtered_feature_bc_matrix ;

    # create dummy
    touch ${meta.cellranger_prefix}/outs/filtered_feature_bc_matrix/matrix.h5ad.gz
    """
}
