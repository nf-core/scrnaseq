/*
 * Reformat design file and check validity
 */
process GENE_MAP {
    tag "$gtf"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path gtf

    output:
    path "transcripts_to_genes.txt" , emit: gene_map

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
