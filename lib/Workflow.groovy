/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.3901628\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    // Exit pipeline if incorrect --genome key provided
    static void genomeExists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(1)
        }
    }

    /*
     * Get workflow summary for MultiQC
     */
    static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    /*
    * Format the protocol
    * Given the protocol paramter (params.protocol) and the aligner (params.aligner),
    * this function formats the protocol such that it is fit for the respective
    * subworkflow
    */
    static formatProtocol(protocol, aligner) {
        String new_protocol = protocol
        String chemistry = ""
        
        // alevin
        if (aligner == "alevin") {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "chromium"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "chromium"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "chromiumV3"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "dropseq"
            }
        }

        // star
        else if (aligner == "star") {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "CB_UMI_Simple"
                    break
                case "smartseq":
                    new_protocol = "SmartSeq"
            }
        }

        // kallisto bustools
        else {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "10XV1"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "10XV2"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "10XV3"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "DROPSEQ"
                    break
                case "smartseq":
                    new_protocol = "SMARTSEQ"
            }
        }
        else {
        exit 1, "Aligner not recognized."
        }

        return [new_protocol, chemistry]
    }
}
