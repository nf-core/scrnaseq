process SALMON_ALEVIN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::salmon=1.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.4.0--h84f40af_1' :
        'quay.io/biocontainers/salmon:1.4.0--h84f40af_1' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
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

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    salmon alevin \\
        -l ISR \\
        -p $task.cpus \\
        -1 ${forward.join( " " )} \\
        -2 ${reverse.join( " " )} \\
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
