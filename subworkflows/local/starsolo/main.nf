/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { STAR_ALIGN  } from '../../../modules/local/star_align'
include { STAR_GENOMEPARAMS_UPGRADE } from '../../../modules/local/star_genomeparams_upgrade'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { STAR_GENOMEGENERATE }         from '../../../modules/nf-core/star/genomegenerate/main'


workflow STARSOLO {
    take:
    genome_fasta
    gtf
    star_index               // path: /path/to/star/index/ (or null)
    star_index_legacy        // boolean: upgrade STAR 2.6.x genomeParameters.txt to 2.7.4a schema
    protocol
    barcode_whitelist
    ch_fastq
    star_feature
    other_10x_parameters

    main:
    ch_versions = channel.empty()

    assert star_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

    assert gtf: "Must provide a gtf file ('--gtf') for STARSOLO"

    /*
    * Build STAR index if not supplied, or upgrade legacy iGenomes metadata when requested
    */
    if (!star_index) {
        STAR_GENOMEGENERATE(
            genome_fasta.map{ f -> [[id: f.baseName], f]},
            gtf.map{ g -> [[id: g.baseName], g]}
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index.collect()
    }
    else {
        // Pre-built STAR index supplied by the user. When star_index_legacy is set
        // (genomes-map opt-in for indices built with STAR 2.6.x, e.g. AWS iGenomes),
        // route through STAR_GENOMEPARAMS_UPGRADE to rewrite `versionGenome 20201` and
        // add the genomeType / genomeTransformType / genomeTransformVCF fields that
        // STAR 2.7.4a+ requires. Modern indices skip the adapter entirely.
        def ch_star_raw = channel.value([ [:], file(star_index, checkIfExists: true) ])
        if (star_index_legacy) {
            log.warn(
                "Using a legacy AWS iGenomes STAR index. nf-core/scrnaseq will update the STAR metadata for " +
                "compatibility, but we recommend regenerating the index for production runs. See " +
                "https://nf-co.re/scrnaseq/dev/docs/usage#reference-genome-options"
            )
            STAR_GENOMEPARAMS_UPGRADE(ch_star_raw)
            ch_star_index = STAR_GENOMEPARAMS_UPGRADE.out.index
        }
        else {
            ch_star_index = ch_star_raw
        }
    }

    /*
    * Perform mapping with STAR
    */
    STAR_ALIGN(
        ch_fastq,
        ch_star_index,
        gtf,
        barcode_whitelist,
        protocol,
        star_feature,
        other_10x_parameters
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    raw_counts = STAR_ALIGN.out.raw_counts
        .join(STAR_ALIGN.out.raw_velocyto, remainder: true)
        .map{
            meta, count, velocity ->
                [meta + [input_type: 'raw'], velocity ? [count, velocity] : [count]]
        }

    filtered_counts = STAR_ALIGN.out.filtered_counts
        .join(STAR_ALIGN.out.filtered_velocyto, remainder: true)
        .map{ meta, count, velocity ->
            [meta + [input_type: 'filtered'], velocity ? [count, velocity] : [count]]
        }

    emit:
    ch_versions
    // get rid of meta for star index
    star_result     = STAR_ALIGN.out.tab
    star_counts     = STAR_ALIGN.out.counts
    raw_counts      = raw_counts
    filtered_counts = filtered_counts
    for_multiqc     = STAR_ALIGN.out.log_final.map{ _meta, logFinal -> logFinal }
}
