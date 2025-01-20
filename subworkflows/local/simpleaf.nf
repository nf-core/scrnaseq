/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { ALEVINQC              } from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        } from '../../../modules/modules/nf-core/simpleaf/index'
include { SIMPLEAF_QUANT        } from '../../../modules/modules/nf-core/simpleaf/quant'

def multiqc_report    = []

workflow SCRNASEQ_SIMPLEAF {

    take:
    genome_fasta
    gtf
    transcript_fasta
    simpleaf_index
    txp2gene
    barcode_whitelist
    resolution
    ch_fastq
    map_dir

    main:
    ch_versions = Channel.empty()

    // we have four types of input:
    // 1. genome fasta and gtf -> build augmented index including spliced transcripts and intronic information
    // 2. transcript fasta and t2g file -> build index directly from transcript fasta
    // 3. simpleaf index and t2g file -> use the index directly
    // 4. mapping results -> skip mapping
    assert ( txp2gene && simpleaf_index  && !transcript_fasta && !genome_fasta && !gtf&& !map_dir ) || // existing index
           ( txp2gene && !simpleaf_index && transcript_fasta && !genome_fasta && !gtf && !map_dir ) || // transcript fasta
           ( !simpleaf_index && !transcript_fasta && genome_fasta && gtf && !map_dir ) || // genome fasta and gtf
           ( !simpleaf_index && !transcript_fasta && !genome_fasta && !gtf && map_dir ) : // existing mapping
        """Must provide either of the followings: (i) --genome_fasta + --gtf, (ii) --transcript_fasta + --txp2gene, (iii) --simpleaf_index + --txp2gene, and (iv) --map_dir""".stripIndent()

    /*
    * Build salmon index
    */
    if ( !simpleaf_index && !map_dir ) {
        SIMPLEAF_INDEX( genome_fasta, gtf, transcript_fasta )
        simpleaf_index = SIMPLEAF_INDEX.out.index.collect()
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
        simpleaf_index,
        txp2gene,
        resolution,
        barcode_whitelist,
        map_dir // this takes a directory of mapping result. Not applicable here
    )
    ch_versions = ch_versions.mix(SIMPLEAF_QUANT.out.versions)

    /*
    * Run alevinQC
    */
    ALEVINQC( SIMPLEAF_QUANT.out.quant, SIMPLEAF_QUANT.out.quant, SIMPLEAF_QUANT.out.map )
    ch_versions = ch_versions.mix(ALEVINQC.out.versions)

    emit:
    ch_versions
    simpleaf_results = SIMPLEAF_QUANT.out.simpleaf.map{ meta, files -> [meta + [input_type: 'raw'], files] }
    alevinqc       = ALEVINQC.out.report
}
