process ANNDATAR_CONVERT {
    tag "${meta.id}"

    label 'process_medium'

    container "fmalmeida/anndatar:dev"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}_standardized.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.id}_standardized.Rds"), emit: rds

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'anndatar_convert.R'

    stub:
    """
    touch ${meta.id}_standardized.h5ad
    touch ${meta.id}_standardized.Rds
    """
}
