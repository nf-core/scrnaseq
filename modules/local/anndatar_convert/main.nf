process ANNDATAR_CONVERT {

    //
    // This module uses the anndata R package to convert h5ad files in different formats
    //

    tag "${meta.id}"

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d8036a262b63ab783572e3cb3120b6b138593e19d93a1059b9bfef80785c7218/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_bioconductor-singlecellexperiment_r-seurat:a0f51df063bb9b2a' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}_${meta.input_type}_matrix*.rds"), emit: rds
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'anndatar_convert.R'

    stub:
    """
    touch ${meta.id}_${meta.input_type}_matrix.Rds
    touch versions.yml
    """
}
