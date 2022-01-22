/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { MULTIQC }                     from '../../modules/local/multiqc_alevin'
include { STAR_ALIGN }                  from '../../modules/local/star_align'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/modules/star/genomegenerate/main'


def multiqc_report    = []

workflow STARSOLO {
    take:
    genome_fasta
    gtf
    star_index
    protocol
    barcode_whitelist
    ch_fastq

    main:
    ch_versions = Channel.empty()

    assert star_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf') if no index is given!"

    assert gtf: "Must provide a gtf file ('--gtf') for STARSOLO"

    /*
    * Build STAR index if not supplied
    */
    if (!params.star_index) {
        STAR_GENOMEGENERATE( genome_fasta, gtf )
        star_index = STAR_GENOMEGENERATE.out.index.collect()
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    /*
    * Perform mapping with STAR
    */
    STAR_ALIGN(
        ch_fastq,
        star_index,
        gtf,
        barcode_whitelist,
        protocol
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
    ch_star_multiqc = STAR_ALIGN.out.log_final


    // /*
    // * MultiQC
    // */
    // if (!params.skip_multiqc) {
    //     workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //         ch_star_multiqc.collect{it[1]}.ifEmpty([]),
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }

    emit:
    ch_versions
    star_result = STAR_ALIGN.out.tab

}
