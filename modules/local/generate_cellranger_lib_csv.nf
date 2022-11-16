process GENERATE_LIB_CSV {
    tag "$samplesheet"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(multi_meta), path(fastqs)

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/scrnaseq/bin/
    def multi_meta_info = multi_meta.collate(2).transpose()
    def sample_types = multi_meta_info[0].join(",")
    def sample_names = multi_meta_info[1].join(",")
    def lib_csv = meta.id + "_lib.csv"
    """
    generate_lib_csv.py \\
        --sample_types $sample_types \\
        --sample_names $sample_names \\
        --out $lib_csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
