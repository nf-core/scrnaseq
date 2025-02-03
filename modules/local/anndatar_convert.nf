process ANNDATAR_CONVERT {

    //
    // This module uses the anndata R package to convert h5ad files in different formats
    //

    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/nfcore/anndatar:20241129'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}_${meta.input_type}_matrix*.rds"), emit: rds
    path  "versions.yml",                                               emit: versions

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
