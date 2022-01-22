/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GENE_MAP }                          from '../../modules/local/gene_map'
include { KALLISTOBUSTOOLS_COUNT }            from '../../modules/local/kallistobustools_count'
include { MULTIQC }                           from '../../modules/local/multiqc_kb'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'
include { KALLISTOBUSTOOLS_REF }        from '../../modules/nf-core/modules/kallistobustools/ref/main'

def multiqc_report    = []

workflow KALLISTO_BUSTOOLS {
    take:
    genome_fasta
    gtf
    kallisto_index
    txp2gene
    protocol
    chemistry
    kb_workflow
    ch_fastq

    main:
    ch_versions = Channel.empty()

    assert kallisto_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf') if no index is given!"

    assert txp2gene || gtf:
        "Must provide either a GTF file ('--gtf') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"

    /*
    * Generate Kallisto Gene Map if not supplied and index is given
    * If no index is given, the gene map will be generated in the 'kb ref' step
    */
    if (!txp2gene && kallisto_index) {
        GENE_MAP( gtf )
        txp2gene = GENE_MAP.out.gene_map
        ch_versions = ch_versions.mix(GENE_MAP.out.versions)
    }

    /*
    * Generate kallisto index
    */
    if (!kallisto_index) {
        KALLISTOBUSTOOLS_REF( genome_fasta, gtf, kb_workflow )
        txp2gene = KALLISTOBUSTOOLS_REF.out.t2g.collect()
        kallisto_index = KALLISTOBUSTOOLS_REF.out.index.collect()
        ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_REF.out.versions)
    }

    /*
    * Quantification with kallistobustools count
    */
    KALLISTOBUSTOOLS_COUNT(
        ch_fastq,
        kallisto_index,
        txp2gene,
        [],
        [],
        false,
        false,
        kb_workflow,
        protocol
    )
    ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_COUNT.out.versions)



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
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }

    emit:
    ch_versions
    counts = KALLISTOBUSTOOLS_COUNT.out.counts


}
