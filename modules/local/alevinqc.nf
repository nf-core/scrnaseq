process ALEVINQC {

    //
    // This module executes alevinfry QC reporting tool on alevin-fry results
    //

    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-alevinqc=1.18.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-alevinqc:1.18.0--r43hf17093f_0' :
        'biocontainers/bioconductor-alevinqc:1.18.0--r43hf17093f_0' }"

    // all metas are the same
    input:
    tuple val(meta), path(quant_dir)
    tuple val(meta1), path(permit_dir)
    tuple val(meta2), path(map_dir)

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
        mapDir = "${map_dir}",
        permitDir= "${permit_dir}",
        quantDir = "${quant_dir}",
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
