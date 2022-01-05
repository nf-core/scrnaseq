process CELLRANGER_MKREF {

    label 'process_medium'

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
    path gtf
    path fasta

    output:
    path 'reference_genome' , emit: reference_genome
    path '*.version.txt'    , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    cellranger mkref \\
        --genome=reference_genome \\
        --genes=${gtf} \\
        --fasta=${fasta} \\
        --nthreads=${task.cpus}

    echo \$(cellranger --version 2>&1) | sed 's/^.*cellranger //; s/ .*\$//' > ${software}.version.txt
    """
}
