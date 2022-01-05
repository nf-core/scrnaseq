process CELLRANGER_MKGTF {

    label 'process_low'

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
    path gtf

    output:
    path '*.gtf'            , emit: gtf
    path '*.version.txt'    , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    cellranger mkgtf \\
        ${gtf} \\
        ${gtf.baseName}_mkgtf.gtf \\
        ${options.args}

    echo \$(cellranger --version 2>&1) | sed 's/^.*cellranger //; s/ .*\$//' > ${software}.version.txt
    """
}
