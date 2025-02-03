/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { GFFREAD_TRANSCRIPTOME } from '../../modules/local/gffread_transcriptome'
include { ALEVINQC              } from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        } from '../../modules/local/simpleaf_index'
include { SIMPLEAF_QUANT        } from '../../modules/local/simpleaf_quant'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP                      } from '../../modules/nf-core/gunzip/main'
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/gffread/main'

// def multiqc_report    = []

workflow SCRNASEQ_ALEVIN {

    take:
    genome_fasta
    gtf
    transcript_fasta
    salmon_index
    txp2gene
    barcode_whitelist
    protocol
    ch_fastq


    main:
    ch_versions = Channel.empty()

    assert (salmon_index && txp2gene) || (genome_fasta && gtf) || (transcript_fasta && txp2gene):
        """Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf'),
        or a genome fasta file and a transcriptome fasta file ('--transcript_fasta`)
        if no index and txp2gene is given!"""
        .stripIndent()

    /*
    * Build salmon index
    */
    if (!salmon_index) {
        SIMPLEAF_INDEX( genome_fasta, transcript_fasta, gtf )
        salmon_index = SIMPLEAF_INDEX.out.index.collect()
        transcript_tsv = SIMPLEAF_INDEX.out.transcript_tsv.collect()
        ch_versions = ch_versions.mix(SIMPLEAF_INDEX.out.versions)

        if (!txp2gene) {
            txp2gene = transcript_tsv
        }
    }

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
    alevin_results = SIMPLEAF_QUANT.out.alevin_results.map{ meta, files -> [meta + [input_type: 'raw'], files] }
    alevinqc       = ALEVINQC.out.report
}
