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

    // simple loop to separate forward from reverse pairs
    def forward_pairs = []
    def reverse_pairs = []
    for (n in 1..reads.toList().size()) {
        current_index = n - 1
        if ( n % 2 == 0 ) { reverse_pairs.add(reads[current_index]) }
        else { forward_pairs.add(reads[current_index]) }
    }
    """
    salmon alevin \\
        -l ISR \\
        -p $task.cpus \\
        -1 ${forward_pairs.join( " " )} \\
        -2 ${reverse_pairs.join( " " )} \\
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
