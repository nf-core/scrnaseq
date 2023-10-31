process CELLRANGERARC_GENERATECONFIG {
    tag "$samplesheet"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val(fasta)
    val(gtf)
    val(motifs)

    output:
    path '*config'     , emit: config
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/scrnaseq/bin/
    def args = task.ext.args ?: ''
    """
    generate_config.py \\
        --fasta $fasta \\
        --gtf $gtf \\
        --motifs $motifs \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
