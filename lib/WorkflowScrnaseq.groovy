//
// This file holds several functions specific to the workflow/scrnaseq.nf in the nf-core/scrnaseq pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowScrnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExists(params, log)

        if (!params.input) {
            log.error "Please provide an input samplesheet with --input"
            System.exit(1)
        }

        if (params.aligner == "cellranger-atac" || params.aligner == "cellranger-arc") {
            if(!params.cellranger_index && !params.reference_config){
                log.error "Reference index or config not specified for cellranger-atac or cellranger-arc with e.g. '--cellranger_index index' or '--reference_config config'."
                System.exit(1)
            }           
        }

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += '    </dl>\n'
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += 'data: |\n'
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }//
    // Exit pipeline if incorrect --genome key provided
    static void genomeExists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    /*
    * Format the protocol
    * Given the protocol paramter (params.protocol) and the aligner (params.aligner),
    * this function formats the protocol such that it is fit for the respective
    * subworkflow
    */
    static formatProtocol(protocol, aligner) {
        String new_protocol = protocol
        String chemistry = ''
        String other_parameters = ''

        // alevin
        if (aligner == 'alevin') {
            switch (protocol) {
                case '10XV1':
                    new_protocol = '10xv1'
                    chemistry = 'V1'
                    break
                case '10XV2':
                    new_protocol = '10xv2'
                    chemistry = 'V2'
                    break
                case '10XV3':
                    new_protocol = '10xv3'
                    chemistry = 'V3'
                    break
                // case 'dropseq':
                //     new_protocol = 'dropseq'
            }
        }

        // star
        else if (aligner == 'star') {
            switch (protocol) {
                case '10XV1':
                    new_protocol = 'CB_UMI_Simple'
                    chemistry = 'V1'
                    other_parameters = '--soloUMIlen 10'
                    break
                case '10XV2':
                    new_protocol = 'CB_UMI_Simple'
                    chemistry = 'V2'
                    other_parameters = '--soloUMIlen 10'
                    break
                case '10XV3':
                    new_protocol = 'CB_UMI_Simple'
                    chemistry = 'V3'
                    other_parameters = '--soloUMIlen 12'
                    break
                case 'dropseq':
                    new_protocol = 'CB_UMI_Simple'
                    break
                case 'smartseq':
                    new_protocol = 'SmartSeq'
            }
        }

        // kallisto bustools
        else if (aligner = 'kallisto' ) {
            switch (protocol) {
                case '10XV1':
                    new_protocol = '10XV1'
                    chemistry = 'V1'
                    break
                case '10XV2':
                    new_protocol = '10XV2'
                    chemistry = 'V2'
                    break
                case '10XV3':
                    new_protocol = '10XV3'
                    chemistry = 'V3'
                    break
                case 'dropseq':
                    new_protocol = 'DROPSEQ'
                    break
                case 'smartseq':
                    new_protocol = 'SMARTSEQ'
            }
        }
        else {
            exit 1, 'Aligner not recognized.'
        }

        return [new_protocol, chemistry, other_parameters]
    }

}
