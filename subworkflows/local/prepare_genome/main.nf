/*
 * Prepare reference FASTA and GTF for alignment (gunzip, filter to genome sequences, optional GTF source fix)
 */

include { GUNZIP as GUNZIP_FASTA              } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                } from '../../../modules/nf-core/gunzip/main'
include { CUSTOM_GTFFILTER as GTF_GENE_FILTER } from '../../../modules/nf-core/custom/gtffilter/main'
include { GAWK as GTF_SOURCE_FIX              } from '../../../modules/nf-core/gawk/main'

workflow PREPARE_GENOME {
    take:
    fasta
    gtf
    gtf_source_fix

    main:
    ch_fasta    = []
    ch_gtf      = []

    if (fasta) {
        fasta_file = file(fasta, checkIfExists: true)
        ch_fasta = channel.value([[id: fasta_file.baseName], fasta_file])

        if (fasta.endsWith('.gz')) {
            ch_fasta = GUNZIP_FASTA(ch_fasta).gunzip
        }
    }

    if (gtf) {
        gtf_file = file(gtf, checkIfExists: true)
        ch_gtf = channel.value([[id: gtf_file.baseName], gtf_file])

        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF(ch_gtf).gunzip
        }

        if (fasta) {
            GTF_GENE_FILTER(
                ch_gtf,
                ch_fasta
            )
            ch_gtf = GTF_GENE_FILTER.out.gtf
        }


        if (gtf_source_fix) {
            // iGenomes GTF annotations with spaces in the source column (e.g. NCBI GRCh38
            // "Curated Genomic") fail Cell Ranger 10 mkref. Opt-in per genome via
            // gtf_source_has_spaces in the genomes map; see usage docs.
            log.warn(
                "Using an iGenomes GTF with spaces in the source column. nf-core/scrnaseq will rewrite the GTF source field for " +
                "Cell Ranger compatibility, but we recommend current reference annotations for production runs. See " +
                "https://nf-co.re/scrnaseq/dev/docs/usage#reference-genome-options"
            )

            GTF_SOURCE_FIX(ch_gtf, [], false)
            ch_gtf = GTF_SOURCE_FIX.out.output
        }
    }

    emit:
    fasta      = ch_fasta.collect()
    gtf        = ch_gtf.collect()
}
