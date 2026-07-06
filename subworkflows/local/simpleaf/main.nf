/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { QCATCH                } from '../../../modules/nf-core/qcatch/main'
include { SIMPLEAF_INDEX        } from '../../../modules/nf-core/simpleaf/index'
include { SIMPLEAF_QUANT        } from '../../../modules/nf-core/simpleaf/quant'

workflow SIMPLEAF {

    take:
    ch_genome_fasta // channel
    ch_genome_gtf   // channel
    transcript_fasta
    simpleaf_index
    txp2gene
    barcode_whitelist
    chemistry
    qcatch_chemistry  // chemistry format for qcatch QC
    skip_qcatch       // skip qcatch QC
    resolution
    ch_fastq   // channel
    map_dir

    main:
    ch_versions = channel.empty()

    /*
    * Build simpleaf index if needed
    * If simpleaf_index is provided, we skip this step
    * If map_dir is provided, we skip this step
    * Otherwise, we build the index
    */
    if ( !simpleaf_index ) {
        if ( !map_dir ) {
            // We do not have a simpleaf index or a map dir, so we need to build the index

            // define input channels for index building
            // we can either use a genome fasta and gtf file pair or a transcript fasta file
            if ( transcript_fasta ) {
                ch_genome_fasta_gtf = [ [:],[],[] ] // meta, genome fasta, genome gtf
                ch_transcript_fasta = [ [id: "${transcript_fasta.getName()}"], transcript_fasta ] // meta, transcript fasta
                if ( !txp2gene ) {
                    error "txp2gene file is required when using `transcript_fasta` to build the index"
                }
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
            // Channel of version
            ch_versions = ch_versions.mix( SIMPLEAF_INDEX.out.versions )

            // ensure txp2gene is a Channel
            if (!txp2gene) {
                txp2gene = SIMPLEAF_INDEX.out.t2g.collect().map { _meta, it -> it }
            } else {
                txp2gene = channel.of( txp2gene )
            }
        } else {
            // we have a map dir, so we do not need to build the index
            simpleaf_index = channel.of( [ [:], [] ] )
        }
    } else {
        // we have a simpleaf index, we use it directly
        // ensure simpleaf index and txp2gene are Channels
        simpleaf_index = channel.of( [ [ id: simpleaf_index.getName() ], simpleaf_index ] )

        // channel or null
        if (txp2gene) {
            txp2gene = channel.of( txp2gene )
        }
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
        // as empty channel cannot be combined, we need this if
        if ( txp2gene ) {
            ch_index_t2g = simpleaf_index.combine( txp2gene ).collect()
        } else {
            ch_index_t2g = simpleaf_index.map {meta, dir -> [meta, dir, []]}.collect()
        }
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
    * Run qcatch QC (optional)
    * qcatch_chemistry may be null for protocols without a QCatch chemistry mapping
    * (e.g. 10XV1 and dropseq); the module then omits --chemistry and QCatch infers
    * it from metadata or uses --n_partitions (set via --qcatch_n_partitions).
    */
    ch_qcatch_report = channel.empty()
    if ( !skip_qcatch ) {
        // Map quant channel to include chemistry for qcatch: tuple(meta, chemistry, quant_dir)
        ch_qcatch_input = ch_af_quant.map { meta, quant_dir -> [meta, qcatch_chemistry, quant_dir] }
        QCATCH( ch_qcatch_input )
        // QCATCH emits versions via the `versions` topic (collected in the main workflow)
        ch_qcatch_report = QCATCH.out.report
    }


    emit:
    ch_versions
    txp2gene
    index       = simpleaf_index
    map         = ch_af_map
    quant       = ch_af_quant
    qcatch_report = ch_qcatch_report
}
