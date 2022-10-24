process CELLRANGER_ARC_COUNT {
    tag "$sample"
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger Arc tool. Please use docker or singularity containers."
    }
    container "heylf/cellranger-arc"

    input:
    path library
    path reference

    output:
    tuple val(sample), path("sample-${sample}/outs/*"), emit: outs
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def library_arg = library.name
    def reference_name = reference.name

    """
    cellranger-arc \\
        count \\
        --id='sample-${sample}' \\
        --libraries=$library_arg \\
        --reference=$reference_name \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "sample-${sample}/outs/"
    touch sample-${sample}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
