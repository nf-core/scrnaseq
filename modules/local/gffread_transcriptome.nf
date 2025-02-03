process GFFREAD_TRANSCRIPTOME {

    //
    // This module uses gffread to filter input to generate a transcripts fasta
    //

    tag "${genome_fasta}"
    label 'process_low'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hd03093a_1' :
        'biocontainers/gffread:0.12.7--hd03093a_1' }"

    input:
    path genome_fasta
    path gtf

    output:
    path "${genome_fasta}.transcriptome.fa", emit: transcriptome_extracted
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread -F $gtf -w "${genome_fasta}.transcriptome.fa" -g $genome_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${genome_fasta}.transcriptome.fa
    touch versions.yml
    """
}
