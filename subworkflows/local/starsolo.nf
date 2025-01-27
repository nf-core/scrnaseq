/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { STAR_ALIGN  } from '../../modules/local/star_align'
include { MTX_TO_H5AD } from '../../modules/local/mtx_to_h5ad'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { GUNZIP }                      from '../../modules/nf-core/gunzip/main'
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/star/genomegenerate/main'


def multiqc_report    = []

workflow STARSOLO {
    take:
    genome_fasta
    gtf
    star_index
    protocol
    barcode_whitelist
    ch_fastq
    star_feature
    other_10x_parameters

    main:
    ch_versions = Channel.empty()

    assert star_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

    assert gtf: "Must provide a gtf file ('--gtf') for STARSOLO"

    /*
    * Build STAR index if not supplied
    */
    if (!star_index) {
        STAR_GENOMEGENERATE(
            genome_fasta.map{ f -> [[id: f.baseName], f]},
            gtf.map{ g -> [[id: g.baseName], g]}
        )
        star_index = STAR_GENOMEGENERATE.out.index.collect()
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    /*
    * Perform mapping with STAR
    */
    STAR_ALIGN(
        ch_fastq,
        star_index,
        gtf,
        barcode_whitelist,
        protocol,
        star_feature,
        other_10x_parameters
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    raw_counts = STAR_ALIGN.out.raw_counts
        .join(STAR_ALIGN.out.raw_velocyto, remainder: true)
        .map{ meta, count, velocity -> {
                return [meta + [input_type: 'raw'], velocity ? [count, velocity] : [count]]
            }
        }

    filtered_counts = STAR_ALIGN.out.filtered_counts
        .join(STAR_ALIGN.out.filtered_velocyto, remainder: true)
        .map{ meta, count, velocity -> {
                return [meta + [input_type: 'filtered'], velocity ? [count, velocity] : [count]]
            }
        }

    emit:
    ch_versions
    // get rid of meta for star index
    star_result     = STAR_ALIGN.out.tab
    star_counts     = STAR_ALIGN.out.counts
    raw_counts      = raw_counts
    filtered_counts = filtered_counts
    for_multiqc     = STAR_ALIGN.out.log_final.map{ meta, it -> it }
}
