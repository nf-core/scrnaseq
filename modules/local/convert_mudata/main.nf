process CONVERT_MUDATA  {
    tag "$meta.id"
    label 'process_single'

    container = 'quay.io/biocontainers/scirpy:0.20.1--pyhdfd78af_0'

    
    input:
    tuple val(meta), path(input_h5ad)
    tuple val(meta), path (input_vdj)

    output:
    tuple val(meta), path("*.mudata.h5mu") , emit: h5mu
    path "versions.yml",  emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def ai = input_vdj ? "-ai $input_vdj" : ''
    
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    convert_mudata.py -ad $input_h5ad $ai

    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    convert_mudata.py --version >> versions.yml
    """

    stub:
    """
    touch matrix.mudata.h5mu
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    convert_mudata.py --version >> versions.yml

    """
}
