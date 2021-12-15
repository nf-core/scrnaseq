include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CELLRANGER_COUNT {

    label 'process_high'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
    tuple val(meta), path('fastqs/*')
    path reference_genome

    output:
    tuple val(meta), path("${prefix}_cellranger")   , emit: cellranger_out
    tuple val(meta), path("${prefix}/*")            , emit: read_counts
    path '*.version.txt'                            , emit: version

    script:
    def prefix = meta.run ? "${meta.sample_name}_${meta.run}" : "${meta.sample_name}"
    def software = getSoftwareName(task.process)

    """
    cellranger count \\
        --id='${prefix}_cellranger' \\
        --fastqs='fastqs' \\
        --sample=${meta.sample_id} \\
        --transcriptome=${reference_genome} \\
        ${options.args}

    mkdir ${prefix}
    cp ${prefix}_cellranger/outs/filtered_feature_bc_matrix/* ${prefix}

    echo \$(cellranger --version 2>&1) | sed 's/^.*cellranger //; s/ .*\$//' > ${software}.version.txt
    """
}
