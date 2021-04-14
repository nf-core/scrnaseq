// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process GFFREAD_TRANSCRIPTOME {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gffread:0.12.1--h2e03b76_1"
    } else {
        container "quay.io/biocontainers/gffread:0.12.1--h2e03b76_1"
    }

    input:
    path bam
    path genome_fasta
    path gtf

    output:n
    path "${genome_fasta}.transcriptome.fa" , emit transcriptome_extracted
    path "*.version.txt"                    , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    gffread -F $gtf -w "${genome_fasta}.transcriptome.fa" -g $genome_fasta

    echo \$(gffread --version 2>&1) > ${software}.version.txt
    """
}
