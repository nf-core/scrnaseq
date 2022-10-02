process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    conda (params.enable_conda ? 'bioconda::salmon=1.5.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0' :
        'quay.io/biocontainers/salmon:1.5.2--h84f40af_0' }"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon"       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kmer_argmatch = args =~ /\-k *(\d+)/
    def k = kmer_argmatch ? kmer_argmatch[0][1] : 31
    """
    grep '^>' $genome_fasta \\
    | cut -d ' ' -f 1 \\
    | sed 's/>//g' > decoys.txt

    cat $genome_fasta \\
    | awk '!/^>/ { next } { getline seq } length(seq) >= $k { print \$0 "\\n" seq }' \\
    | gzip -c > gentrome.filtered.fasta.gz

    salmon \\
        index \\
        --threads $task.cpus \\
        -t gentrome.filtered.fasta.gz \\
        -d decoys.txt \\
        $args \\
        -i salmon
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
