/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GFFREAD_TRANSCRIPTOME }             from '../../modules/local/gffread_transcriptome'
include { SALMON_ALEVIN         }             from '../../modules/local/salmon_alevin'
include { ALEVINQC              }             from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        }             from '../../modules/local/simpleaf_index'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/modules/gffread/main'

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

    assert salmon_index || (genome_fasta && gtf) || (genome_fasta && transcript_fasta):
        """Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf'), or a genome fasta file
        and a transcriptome fasta file ('--transcript_fasta`) if no index is given!""".stripIndent()

    assert txp2gene || gtf:
        "Must provide either a GTF file ('--gtf') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"


    /*
    * Build salmon index
    */
    if (!salmon_index) {
        SIMPLEAF_INDEX( genome_fasta, transcript_fasta, gtf )
        salmon_index = SIMPLEAF_INDEX.out.index.collect()
        ch_versions = ch_versions.mix(SIMPLEAF_INDEX.out.versions)
    }

    /*
    * Build txp2gene map
    */
    if (!txp2gene){
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

    /*
    * Run alevinQC
    */
    ALEVINQC( SALMON_ALEVIN.out.alevin_results )
    ch_versions = ch_versions.mix(ALEVINQC.out.versions)

    emit:
    ch_versions
    alevin_results = SALMON_ALEVIN.out.alevin_results
    alevinqc = ALEVINQC.out.report
    for_multiqc = SALMON_ALEVIN.out.alevin_results.collect{it[1]}.ifEmpty([])

}
