process ALEVINQC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bioconductor-alevinqc=1.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-alevinqc:1.10.0--r41hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-alevinqc:1.10.0--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(alevin_results)

    output:
    tuple val(meta), path("alevin_report_${meta.id}.html"), emit: report
    path "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    require(alevinQC)
    alevinQCReport(baseDir = "${alevin_results}", sampleId = "${prefix}",
                outputFile = "alevin_report_${meta.id}.html",
                outputFormat = "html_document",
                outputDir = "./", forceOverwrite = TRUE)
    write(as.character(packageVersion("alevinQC")), paste0("${software}", ".version.txt"))
    """
}
