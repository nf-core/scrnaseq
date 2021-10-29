// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::multiqc=1.10.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.10.1--py_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path ('salmon_alevin/*')

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path "*variants_metrics_mqc.csv", optional:true, emit: csv_variants
    path "*assembly_metrics_mqc.csv", optional:true, emit: csv_assembly
    path "*_plots"                  , optional:true, emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $options.args $custom_config .
    """
}
