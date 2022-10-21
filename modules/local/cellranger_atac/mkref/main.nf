process CELLRANGER_ATAC_MKREF {
    tag "$fasta"
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "heylf/cellranger-atac"

    input:
    val ref_config_name

    output:
    path "${ref_config_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellranger-atac \\
        mkref \\
        --config=$ref_config_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-atac: \$(echo \$( cellranger-atac --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
