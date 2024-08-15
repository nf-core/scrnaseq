process BAMTOFASTQ10X {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/10x_bamtofastq:1.4.1--hdbdd923_2':
        'biocontainers/10x_bamtofastq:1.4.1--hdbdd923_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}/**/*.fastq.gz"), emit: fastq
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamtofastq \\
        $args \\
        $bam \\
        $prefix

    out_dir=\$(find . -type d -maxdepth 2 -print | grep -m1 '${meta.sample_id}_0_1')
    echo \${out_dir}

    for file in $prefix/${meta.sample_id}_0_1*/*.fastq.gz; 
    do
        echo \$file
        mv "\$file" "\${file/bamtofastq/$prefix}"; 
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtofastq10x: \$(bamtofastq --version |& sed '1!d ; s/bamtofastq //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtofastq10x: \$(bamtofastq --version |& sed '1!d ; s/bamtofastq //')
    END_VERSIONS
    """
}
