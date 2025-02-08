/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { ALEVINQC              } from '../../modules/local/alevinqc'
include { SIMPLEAF_INDEX        } from '../../modules/nf-core/simpleaf/index'
include { SIMPLEAF_QUANT        } from '../../modules/nf-core/simpleaf/quant'

workflow SCRNASEQ_SIMPLEAF {

    take:
    ch_genome_fasta // channel
    ch_genome_gtf   // channel
    transcript_fasta
    simpleaf_index
    txp2gene
    barcode_whitelist
    chemistry
    resolution
    ch_fastq   // channel
    map_dir

    main:
    ch_versions = Channel.empty()

    /*
    * Build simpleaf index
    */
    if ( !simpleaf_index || !map_dir ) {
        // define input channels for index building
        // we can either use the genome fasta and gtf files or the transcript fasta file
        if ( transcript_fasta ) {
            ch_genome_fasta_gtf = [ [:],[],[] ] // meta, genome fasta, genome gtf
            ch_transcript_fasta = [ [id: "${transcript_fasta.getName()}"], transcript_fasta ] // meta, transcript fasta
        } else {
            ch_genome_fasta_gtf = ch_genome_fasta.combine( ch_genome_gtf ).map{ fasta, gtf -> [[id: "${fasta.getName()}"], fasta, gtf] }
            ch_transcript_fasta = [ [:], [] ] // meta, transcript fasta
        }

        SIMPLEAF_INDEX(
            ch_genome_fasta_gtf,
            ch_transcript_fasta,
            [[:], []], // meta, probe CSV
            [[:], []] // meta, feature CSV
        )
        // Channel of tuple(meta, index dir)
        simpleaf_index = SIMPLEAF_INDEX.out.index.collect()
        // Channel of t2g path or empty
        t2g = SIMPLEAF_INDEX.out.t2g.collect().map { _meta, it -> it }
        ch_versions = ch_versions.mix( SIMPLEAF_INDEX.out.versions )

        // ensure txp2gene is a Channel
        if (!txp2gene) {
            txp2gene = t2g
        } else {
            txp2gene = Channel.of( txp2gene )
        }
    } else {
        // ensure simpleaf index and txp2gene are Channels
        simpleaf_index = Channel.of( [ [ id: simpleaf_index.getName() ], simpleaf_index ] )
    }

    // define input channels for quantification
    // we can either use the mapping results or the reads and index files
    if ( map_dir ) {
        // meta, chemistry, files
        ch_chemistry_reads = [ [:],[],[] ]
        // meta, index, t2g file
        ch_index_t2g = [ [:],[],[] ]
        // meta, map dir
        ch_map_dir = [ [ id: map_dir.baseName ], map_dir ]
    } else {
        // meta, chemistry, files
        ch_chemistry_reads = ch_fastq.map{ meta, files -> [meta + ["chemistry": chemistry], chemistry, files] }
        // meta, index, t2g file
        ch_index_t2g = simpleaf_index.combine( txp2gene ).collect()
        // meta, map dir
        ch_map_dir = [ [:],[] ]
    }

    /*
    * Perform quantification with simpleaf quant
    */
    SIMPLEAF_QUANT (
        ch_chemistry_reads,
        ch_index_t2g,
        // meta, cell filtering method, cell filtering params, whitelist
        [[:], "unfiltered-pl", [], barcode_whitelist ],
        resolution,
        ch_map_dir
    )
    ch_versions = ch_versions.mix(SIMPLEAF_QUANT.out.versions)

    ch_af_map = map_dir ? ch_map_dir : SIMPLEAF_QUANT.out.map
    ch_af_quant = SIMPLEAF_QUANT.out.quant
    /*
    * Run alevinQC
    */
    ALEVINQC( ch_af_quant, ch_af_quant, ch_af_map )
    ch_versions = ch_versions.mix(ALEVINQC.out.versions)


    emit:
    ch_versions
    txp2gene
    index       = simpleaf_index
    map         = ch_af_map
    quant       = ch_af_quant
    alevinqc    = ALEVINQC.out.report
}
