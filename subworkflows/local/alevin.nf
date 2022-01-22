/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GFFREAD_TRANSCRIPTOME }             from '../../modules/local/gffread_transcriptome'
include { SALMON_ALEVIN }                     from '../../modules/local/salmon_alevin'
include { ALEVINQC }                          from '../../modules/local/alevinqc'
include { MULTIQC }                           from '../../modules/local/multiqc_alevin'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/modules/gffread/main'
include { SALMON_INDEX }                from '../../modules/nf-core/modules/salmon/index/main'

def multiqc_report    = []

workflow SCRNASEQ_ALEVIN {

    take:
    genome_fasta
    gtf
    transcript_fasta
    salmon_index
    txp2gene
    barcode_whitelist
    protocol
    chemistry
    ch_fastq


    main:
    ch_versions = Channel.empty()

    assert salmon_index || (genome_fasta && gtf) || (genome_fasta && transcriptome_fasta):
        """Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf'), or a genome fasta file
        and a transcriptome fasta file ('--transcriptome_fasta`) if no index is given!""".stripIndent()

    assert txp2gene || gtf:
        "Must provide either a GTF file ('--gtf') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"


    /*
    * Build salmon index
    */
    if (!params.salmon_index) {
        // Preprocessing - Extract transcriptome fasta from genome fasta
        if (!params.transcript_fasta) {
            GFFREAD_TRANSCRIPTOME( genome_fasta, gtf )
            transcriptome_fasta = GFFREAD_TRANSCRIPTOME.out.transcriptome_extracted
            ch_versions = ch_versions.mix(GFFREAD_TRANSCRIPTOME.out.versions)
        }
        SALMON_INDEX( genome_fasta, transcriptome_fasta )
        salmon_index = SALMON_INDEX.out.index.collect()
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    }

    /*
    * Build txp2gene map
    */
    if (!params.txp2gene){
        GFFREAD_TXP2GENE( gtf )
        txp2gene = GFFREAD_TXP2GENE.out.gtf
        // Only collect version if not already done for gffread
        ch_versions = ch_versions.mix(GFFREAD_TXP2GENE.out.versions)
    }

    /*
    * Perform quantification with salmon alevin
    */
    SALMON_ALEVIN (
        ch_fastq,
        salmon_index,
        txp2gene,
        protocol,
        barcode_whitelist
    )
    ch_versions = ch_versions.mix(SALMON_ALEVIN.out.versions)
    ch_salmon_multiqc = SALMON_ALEVIN.out.alevin_results

    /*
    * Run alevinQC
    */
    ALEVINQC( SALMON_ALEVIN.out.alevin_results )
    ch_versions = ch_versions.mix(ALEVINQC.out.versions)


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
    //         ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }

    emit:
    ch_versions
    alevin_results = SALMON_ALEVIN.out.alevin_results
    alevinqc = ALEVINQC.out.report

}
