process SALMON_ALEVIN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::salmon=1.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.4.0--h84f40af_1' :
        'quay.io/biocontainers/salmon:1.4.0--h84f40af_1' }"

    input:
    tuple val(meta), path(reads)
    path index
    path txp2gene
    val protocol
    path whitelist

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
    path  "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    salmon alevin \\
        -l ISR \\
        -p $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --${protocol} \\
        -i $index \\
        --tgMap $txp2gene \\
        --dumpFeatures --dumpMtx \\
        $args \\
        -o ${prefix}_alevin_results

    gzip -cdf ${whitelist} > ${prefix}_alevin_results/alevin/whitelist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
