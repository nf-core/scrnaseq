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
        tuple val(meta), path("${prefix}_cellranger"), emit: cellranger_out
        tuple val(meta), path("${prefix}/*"), emit: read_counts

    script:
        prefix = meta.run ? "${meta.sample_name}_${meta.run}" : "${meta.sample_name}"

        cellranger_count_command = "cellranger count --id='${prefix}_cellranger' --fastqs='fastqs' --sample=${meta.sample_id} --transcriptome=${reference_genome} ${options.args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] cellranger count command: " + cellranger_count_command)
        }

       //SHELL
        """
        ${cellranger_count_command}
        mkdir ${prefix}
        cp ${prefix}_cellranger/outs/filtered_feature_bc_matrix/* ${prefix}
        """
}