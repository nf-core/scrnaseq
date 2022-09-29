process SIMPLEAF_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::simpleaf=0.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.5.1--h9f5acd7_0' :
        'quay.io/biocontainers/simpleaf:0.5.1--h9f5acd7_0' }"

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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    gzip -dcf $whitelist > whitelist.txt
    simpleaf quant \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -i ${index} \\
        -o ${prefix}_alevin_results \\
        -m $txp2gene \\
        -t $task.cpus \\
        -c $protocol \\
        -u whitelist.txt \\
        $args
    
    mv whitelist.txt ${prefix}_alevin_results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
