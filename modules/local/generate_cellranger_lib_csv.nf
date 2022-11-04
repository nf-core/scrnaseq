process GENERATE_LIB_CSV {
    tag "$samplesheet"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(folders)

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/scrnaseq/bin/
    def folder_GEX = meta.folders[0]
    def folder_ATAC = meta.folders[1]
    def sample_arg = meta.id
    def lib_csv = sample_arg + "_lib.csv"
    """
    generate_lib_csv.py \\
        --in_GEX $folder_GEX \\
        --in_ATAC $folder_ATAC \\
        --sample $sample_arg \\
        --out $lib_csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
