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
    path  "versions.yml"                      , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    require(alevinQC)
    alevinQCReport(baseDir = "${alevin_results}", sampleId = "${prefix}",
                outputFile = "alevin_report_${meta.id}.html",
                outputFormat = "html_document",
                outputDir = "./", forceOverwrite = TRUE)

    yaml::write_yaml(
        list(
            '${task.process}'=list(
                'alevinqc' = paste(packageVersion('alevinQC'), collapse='.')
            )
        ),
        "versions.yml"
    )
    """
}
