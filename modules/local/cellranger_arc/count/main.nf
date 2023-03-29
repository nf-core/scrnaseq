process CELLRANGER_ARC_COUNT {
    tag "$meta.id"
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "heylf/cellranger-arc:2.0.2"

    input:
    tuple val(meta), val(multi_meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs/*"), emit: outs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_name = reference.name

    def multi_meta_info = multi_meta.collate(2).transpose()
    def sample_types = multi_meta_info[0].join(",")
    def sample_names = multi_meta_info[1].join(",")
    def lib_csv = meta.id + "_lib.csv"

    """
    # The following ugly three commands (mkdir, mv, generate_lib_csv)
    # are required because cellranger-arc only deals with abolsute paths
    if [ ! -d "fastqs" ]; then
        mkdir fastqs
    fi

    mv *.fastq.gz fastqs/

    generate_lib_csv.py \\
        --sample_types $sample_types \\
        --sample_names $sample_names \\
        --fastq_folder \$(readlink -f fastqs)\\
        --out $lib_csv

    cellranger-arc \\
        count \\
        --id='${meta.id}' \\
        --libraries=$lib_csv \\
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
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
