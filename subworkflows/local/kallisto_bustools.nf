/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GENE_MAP }                          from '../../modules/local/gene_map'
include {KALLISTOBUSTOOLS_COUNT }             from '../../modules/nf-core/kallistobustools/count/main'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/gunzip/main'
include { KALLISTOBUSTOOLS_REF }        from '../../modules/nf-core/kallistobustools/ref/main'

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
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

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
        t1c = KALLISTOBUSTOOLS_REF.out.cdna_t2c.ifEmpty{ [] }
        t2c = KALLISTOBUSTOOLS_REF.out.intron_t2c.ifEmpty{ [] }
    }

    /*
    * Quantification with kallistobustools count
    */
    KALLISTOBUSTOOLS_COUNT(
        ch_fastq,
        kallisto_index,
        txp2gene,
        t1c,
        t2c,
        protocol
    )

    ch_versions = ch_versions.mix(KALLISTOBUSTOOLS_COUNT.out.versions)

    emit:
    ch_versions
    counts = KALLISTOBUSTOOLS_COUNT.out.count
    txp2gene


}
