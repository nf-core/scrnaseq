/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GFFREAD_TRANSCRIPTOME }             from '../../modules/local/gffread_transcriptome'
include { ALEVINQC              }             from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        }             from '../../modules/local/simpleaf_index'
include { SIMPLEAF_QUANT        }             from '../../modules/local/simpleaf_quant'

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
        transcript_tsv = SIMPLEAF_INDEX.out.transcript_tsv.collect()
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
    SIMPLEAF_QUANT (
        ch_fastq,
        salmon_index,
        transcript_tsv,
        txp2gene,
        protocol,
        barcode_whitelist
    )
    ch_versions = ch_versions.mix(SIMPLEAF_QUANT.out.versions)

    /*
    * Run alevinQC
    */
    ALEVINQC( SIMPLEAF_QUANT.out.alevin_results )
    ch_versions = ch_versions.mix(ALEVINQC.out.versions)

    emit:
    ch_versions
    alevin_results = SIMPLEAF_QUANT.out.alevin_results
    alevinqc = ALEVINQC.out.report
    for_multiqc = SIMPLEAF_QUANT.out.alevin_results.collect{it[1]}.ifEmpty([])

}
