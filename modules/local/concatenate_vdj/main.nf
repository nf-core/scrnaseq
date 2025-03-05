process CONCATENATE_VDJ   {
    tag "$meta.id"
    label 'process_single'

    container = 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'
    
    input:
    tuple val(meta), path(input_vdj, stageAs: '?/*')

    output:
    tuple val(meta), path("*.vdj.h5ad") , emit: h5ad, optional: true
    path "versions.yml",  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:   
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    concatenate_vdj.py -ai ${input_vdj.join(' ')} -id ${meta.collect{ it.id }.join(' ')}
    
    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    concatenate_vdj.py --version >> versions.yml
    """

    stub:
    """
    touch combined_vdj.h5ad
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    concatenate_vdj.py --version >> versions.yml
    
    """   
}