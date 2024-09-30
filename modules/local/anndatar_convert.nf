process ANNDATAR_CONVERT {
    tag "${meta.id}"

    label 'process_medium'

    container "docker://fmalmeida/anndatar:dev" // TODO: Fix

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_${meta.input_type}_matrix.Rds"), emit: rds

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'anndatar_convert.R'

    stub:
    """
    touch ${meta.id}.Rds
    """
}
