process SIMPLEAF_INDEX {
    tag "$transcript_gtf"
    label "process_medium"

    conda (params.enable_conda ? 'bioconda::simpleaf=0.5.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.5.3--h9f5acd7_0' :
        'quay.io/biocontainers/simpleaf:0.5.3--h9f5acd7_0' }"

    input:
    path genome_fasta
    path transcript_fasta
    path transcript_gtf

    output:
    path "salmon/index"              , emit: index
    path "salmon/ref/*_t2g_3col.tsv" , emit: transcript_tsv
    path "versions.yml"              , emit: versions
    path "salmon"                    , emit: salmon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_inputs = (params.transcript_fasta) ? "--refseq $transcript_fasta" : "--gtf $transcript_gtf"
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf index
    simpleaf \\
        index \\
        --threads $task.cpus \\
        --fasta $genome_fasta \\
        $seq_inputs \\
        $args \\
        -o salmon

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
