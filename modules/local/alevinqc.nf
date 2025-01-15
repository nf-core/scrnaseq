process ALEVINQC {

    //
    // This module executes alevinfry QC reporting tool on alevin results
    //

    tag "$meta.id"
    label 'process_low'

    //The alevinqc 1.14.0 container is broken, missing some libraries - thus reverting this to previous 1.12.1 version
    conda "bioconda::bioconductor-alevinqc=1.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-alevinqc:1.12.1--r41h9f5acd7_0' :
        'biocontainers/bioconductor-alevinqc:1.12.1--r41h9f5acd7_0' }"

    input:
    tuple val(meta), path(alevin_results)

    output:
    tuple val(meta), path("alevin_report_${meta.id}.html"), emit: report
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    require(alevinQC)
    alevinFryQCReport(
        mapDir = "${alevin_results}/af_map",
        quantDir = "${alevin_results}/af_quant",
        permitDir= "${alevin_results}/af_quant",
        sampleId = "${prefix}",
        outputFile = "alevin_report_${meta.id}.html",
        outputFormat = "html_document",
        outputDir = "./",
        forceOverwrite = TRUE
    )

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
