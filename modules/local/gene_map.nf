/*
 * Reformat design file and check validity
 */
process GENE_MAP {
    tag "$gtf"
    label 'process_low'


    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path gtf

    output:
    path "transcripts_to_genes.txt" , emit: gene_map

    when:
    task.ext.when == null || task.ext.when

    script:
    if("${gtf}".endsWith('.gz')){
        name = "${gtf.baseName}"
        unzip = "gunzip -f ${gtf}"
    } else {
        unzip = ""
        name = "${gtf}"
    }
    """
    $unzip
    cat $name | t2g.py --use_version > transcripts_to_genes.txt
    """
}
