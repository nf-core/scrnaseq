/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GFFREAD_TRANSCRIPTOME }             from '../../modules/local/gffread_transcriptome'
include { ALEVINQC              }             from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        }             from '../../modules/local/simpleaf_index'
include { SIMPLEAF_QUANT        }             from '../../modules/local/simpleaf_quant'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/gunzip/main'
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/gffread/main'

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
        """Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf'), or a genome fasta file
        and a transcriptome fasta file ('--transcript_fasta`) if no index is given!""".stripIndent()

    if (transcript_fasta) {
        assert txp2gene:
            "Since a built transcript was provided ('--transcript_fasta'), must also provide a simpleaf gene map ('--txp2gene') to use with simpleaf quant!"
    }


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
    * Select txp2gene map
    */
    if (!txp2gene) { txp2gene = SIMPLEAF_INDEX.out.transcript_tsv }

    /*
    * Perform quantification with salmon alevin
    */
    SIMPLEAF_QUANT (
        ch_fastq,
        salmon_index,
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
