//
// Subworkflow with functionality specific to the nf-core/scrnaseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/scrnaseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/scrnaseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    if (params.aligner == 'cellrangermulti') { // the cellrangermulti sub-workflow logic needs that channels have reads separated by feature_type. Cannot merge all.
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta.feature_type, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta.feature_type, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple( by: [0,1] )
            .map{ id, _type, meta, reads -> [ id, meta, reads ] }
            .map { sheet_row ->
                validateInputSamplesheet(sheet_row)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
    } else if (params.aligner == 'cellrangerarc') { // the cellrangerarc sub-workflow logic needs that channels have a meta, type, subsample, fastqs structure.
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2 || (meta.sample_type == "atac" && !meta.fastq_barcode)) {
                    error("Please check input samplesheet -> cellrangerarc requires both paired-end reads and barcode fastq files: ${meta.id}")
                }
                if (meta.sample_type == "atac") {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2, file(meta.fastq_barcode, checkIfExists: true) ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
            }
            .groupTuple()
            .map { structure_input ->
                cellrangerarcStructure(structure_input)
            }
            .set { ch_samplesheet }
    } else {
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { sheet_row ->
                validateInputSamplesheet(sheet_row)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    // Validate cellranger_multi_barcodes if provided and aligner is cellrangermulti
    if (params.aligner == 'cellrangermulti' && params.cellranger_multi_barcodes) {
        validateCellrangerMultiBarcodes()
    }
}

//
// Validate cellranger_multi_barcodes samplesheet for uniqueness and conditional requirements
//
def validateCellrangerMultiBarcodes() {
    def cellranger_multi_barcodes = file(params.cellranger_multi_barcodes).splitCsv(header: true)

    // Get unique samples from input samplesheet for cross-validation
    def inputSamples = file(params.input).splitCsv(header: true).collect { row -> row.sample }.toSet()

    // Check that at least one barcode column is provided for each row
    // and that each sample uses only one type of barcode
    def rowsWithoutBarcodes = []
    def sampleBarcodeTypes = [:]
    def barcodeSamples = [] as Set
    cellranger_multi_barcodes.eachWithIndex { row, idx ->
        def multiplexed_sample_id = row.multiplexed_sample_id
        def rowNum = idx + 2 // +2 for 1-based indexing and header row

        // Collect unique sample names for cross-validation
        barcodeSamples << row.sample

        def barcodeTypes = []
        if (row.probe_barcode_ids) barcodeTypes << 'probe_barcode_ids'
        if (row.cmo_ids)           barcodeTypes << 'cmo_ids'
        if (row.ocm_ids)           barcodeTypes << 'ocm_ids'

        if (barcodeTypes.isEmpty()) {
            rowsWithoutBarcodes << [row: rowNum, multiplexed_sample_id: multiplexed_sample_id]
        }
        sampleBarcodeTypes[multiplexed_sample_id] = [types: barcodeTypes.toSet(), row: rowNum]
    }

    // Validate that at least one barcode identifier is populated in each row
    if (rowsWithoutBarcodes) {
        def errorDetails = rowsWithoutBarcodes.collect { missing -> "row ${missing.row} (${missing.multiplexed_sample_id})" }.join(', ')
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "The following rows have no barcode identifiers: ${errorDetails}. " +
              "Each row must have exactly one of: 'probe_barcode_ids', 'cmo_ids', or 'ocm_ids'.")
    }

    // Validate that no more than one barcode identifier is populated in each row
    def samplesWithMixedBarcodes = sampleBarcodeTypes.findAll { _multiplexed_sample_id, info -> info.types.size() > 1 }
    if (samplesWithMixedBarcodes) {
        def errorMsg = samplesWithMixedBarcodes.collect { multiplexed_sample_id, info ->
            "'${multiplexed_sample_id}' (row ${info.row}) uses multiple barcode types: ${info.types.join(', ')}"
        }.join('; ')
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "Each multiplexed_sample_id should use only one type of barcode identifier. ${errorMsg}")
    }

    // Validate that samples in cellranger_multi_barcodes exist in the input samplesheet
    def unknownSamples = barcodeSamples - inputSamples
    if (unknownSamples) {
        error("Please check cellranger_multi_barcodes samplesheet -> " +
              "The following sample(s) do not exist in the input samplesheet: ${unknownSamples.join(', ')}. " +
              "The 'sample' column in cellranger_multi_barcodes must match 'sample' values in the input samplesheet.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// cellrangerarc structure for samplesheet channel
//
def cellrangerarcStructure(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Validate that the property "sample_type" is present and has valid values
    def valid_sample_types = ["gex", "atac"]
    def sample_type_ok = metas.collect { meta -> meta.sample_type }.unique().every { st -> st in valid_sample_types }
    if (!sample_type_ok) {
        error("Please check input samplesheet -> The property 'sample_type' is required and can only be 'gex' or 'atac'.")
    }

    // Define a new common meta for all the fastqs in this channel instance
    def sampleMeta = metas[0].clone()
    sampleMeta.remove("sample_type")
    sampleMeta.remove("feature_type")

    // Create a list with all the entries of meta.sample_type
    def sampletypes = metas.collect { meta -> meta.sample_type }

    // Create a list with all the base name of the fastq files
    def subsamples = fastqs.collect { fastq ->
        def match = (fastq[0].baseName =~ /^(.*?)_S\d+_L\d+_R\d+_\d+\.fastq(\.gz)?$/)
        if (!match) {
            error("Filename does not follow the expected FASTQ filename convention (SampleName_S1_L001_R1_001.fastq.gz): ${fastq[0]}")
        }
        return match[0][1]
    }

    return [ sampleMeta, sampletypes, subsamples, fastqs.flatten() ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        } else {
            return null
        }
    } else {
        return null
    }
}

//
// iGenomes GTF annotations with spaces in the GTF source column (e.g. NCBI RefSeq "Curated Genomic")
// are incompatible with Cell Ranger 10 mkref; opt-in per genome via gtf_source_has_spaces.
//
def gtfSourceFixNeeded(aligner, genome, genomes, gtf) {
    def genome_entry = genomes && genome ? genomes[genome] : null
    def cellranger_aligner = aligner in ['cellranger', 'cellrangerarc', 'cellrangermulti']
    def gtf_flagged = genome_entry?.gtf_source_has_spaces as Boolean
    def gtf_from_genome = gtf == genome_entry?.gtf
    return cellranger_aligner && gtf_flagged && gtf_from_genome
}

//
// Decide whether the supplied STAR index needs to be routed through the
// STAR_GENOMEPARAMS_UPGRADE adapter. Fires when the active genomes-map entry
// has `star_legacy = true` (set on every iGenomes entry that ships a
// `star` directory) and the user has not overridden the resolved index with
// their own --star_index.
//
def isStarIndexLegacy(genome, genomes, star_index) {
    def genome_entry = genomes && genome ? genomes[genome] : null
    def star_legacy = genome_entry?.star_legacy as Boolean
    def index_from_genome = star_index == genome_entry?.star
    return star_legacy && index_from_genome
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
