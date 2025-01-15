process PARSE_CELLRANGERMULTI_SAMPLESHEET {

    //
    // This module contains a custom script for checking special cellranger multi samplesheet
    //

    label 'process_low'
    publishDir = [ enabled: false ]

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path samplesheet

    output:
    path "cmo_files/*" , emit: cmo,  optional: true
    path "frna_files/*", emit: frna, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_cellrangermulti.py $samplesheet
    """

    stub:
    """
    mkdir -p cmo_files frna_files
    touch frna_files/test.csv
    touch cmd_files/test.csv
    """
}
