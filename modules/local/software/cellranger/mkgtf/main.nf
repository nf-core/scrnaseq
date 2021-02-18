include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CELLRANGER_MKGTF {

    label 'process_low'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
        path(gtf)

    output:
        path("*.gtf")

    script:
        
        mkgtf_command = "cellranger mkgtf ${gtf} ${gtf.baseName}_mkgtf.gtf ${options.args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] filter gtf command: " + mkgtf_command)
        }

       //SHELL
        """
        ${mkgtf_command}
        """
}