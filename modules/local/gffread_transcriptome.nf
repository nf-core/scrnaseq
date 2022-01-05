process GFFREAD_TRANSCRIPTOME {
    tag "${genome_fasta}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gffread:0.12.1--h2e03b76_1"
    } else {
        container "quay.io/biocontainers/gffread:0.12.1--h2e03b76_1"
    }

    input:
    path genome_fasta
    path gtf

    output:
    path "${genome_fasta}.transcriptome.fa" , emit: transcriptome_extracted
    path "*.version.txt"                    , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    gffread -F $gtf -w "${genome_fasta}.transcriptome.fa" -g $genome_fasta

    echo \$(gffread --version 2>&1) > ${software}.version.txt
    """
}
