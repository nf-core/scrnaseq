process SIMPLEAF_QUANT {

    //
    // This module executes simpleaf to perform quantification with alevin
    //

    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::simpleaf=0.10.0-1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.10.0--h9f5acd7_1' :
        'biocontainers/simpleaf:0.10.0--h9f5acd7_1' }"

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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args_list = args.tokenize()
    def prefix    = task.ext.prefix ?: "${meta.id}"

    //
    // check if users are using one of the mutually excludable parameters:
    //    e.g -k,--knee | -e,--expect-cells | -f, --forced-cells
    //
    unzip_whitelist = ""
    unfiltered_command = ""
    save_whitelist     = ""
    if (!(args_list.any { it in ['-k', '--knee', '-e', '--expect-cells', '-f', '--forced-cells']} || meta.expected_cells)) {
        if (whitelist) {
            if (whitelist.name.endsWith('.gz')) {
                unzip_whitelist = "gzip -dcf $whitelist > whitelist.uncompressed.txt"
            } else {
                unzip_whitelist = "cp $whitelist whitelist.uncompressed.txt"
            }
            unfiltered_command = "-u whitelist.uncompressed.txt"
            save_whitelist     = "mv whitelist.uncompressed.txt ${prefix}_alevin_results/"
        }
    }

    // expected cells
    def expect_cells = meta.expected_cells ? "--expect-cells $meta.expected_cells" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    # export required var
    export ALEVIN_FRY_HOME=.
    export NUMBA_CACHE_DIR=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    $unzip_whitelist
    simpleaf quant \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -i ${index} \\
        -o ${prefix}_alevin_results \\
        -m $txp2gene \\
        -t $task.cpus \\
        -c "$protocol" \\
        $expect_cells \\
        $unfiltered_command \\
        $args

    $save_whitelist
    [[ ! -f ${prefix}_alevin_results/af_quant/all_freq.bin ]] && cp ${prefix}_alevin_results/af_quant/permit_freq.bin ${prefix}_alevin_results/af_quant/all_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
