process EMPTYDROPS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy==1.7.2 conda-forge::python-igraph==0.10.4 conda-forge::leidenalg==0.9.1 bioconda::emptydrops==0.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.mtx"), emit: mtx

    script:
    // TODO:
    // check if inputs meet expectations (input directory)
    // add output directory path
    """
    # filter mtx files
    filter_by_emptydrops.py \\
        --mtx_dir $inputs \\
        --out_dir . \\
        --version ${meta.version} ;
    """

    stub:
    """
    touch ${meta.id}_matrix.h5ad
    """
}
