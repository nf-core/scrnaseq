//
// Include modules
//
include { CELLRANGER_MKGTF } from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF } from "../../modules/nf-core/cellranger/mkref/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI {
    take:
        ch_fasta
        ch_gtf
        ch_fastq
        cellranger_gex_index

    main:
        ch_versions    = Channel.empty()
        def empty_file = file("$projectDir/assets/EMPTY", checkIfExists: true)

        // since we merged all data as a meta, now we have a channel per sample, which
        // every item is a meta map for each data-type
        // now we can split it back for passing as input to the module
        ch_fastq
        .flatten()
        .map{ meta ->
            def meta_clone = meta.clone()
            def data_dict  = meta_clone.find{ it.key == "${meta_clone.feature_type}" }
            fastqs = data_dict?.value
            meta_clone.remove( data_dict?.key )
            [ meta_clone, fastqs ]
        }
        .branch {
            meta, fastq ->
                gex: meta.feature_type == "gex"
                    return [ meta, fastq ]
                vdj: meta.feature_type == "vdj"
                    return [ meta, fastq ]
                ab: meta.feature_type == "ab"
                    return [ meta, fastq ]
                beam: meta.feature_type == "beam"
                    return [ meta, fastq ]
                crispr: meta.feature_type == "crispr"
                    return [ meta, fastq ]
                cmo: meta.feature_type == "cmo"
                    return [ meta, fastq ]
        }
        .set { ch_grouped_fastq }

        //
        // Prepare gex reference (Normal Ref)
        //
        if ( !cellranger_gex_index ) {

            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF ( ch_gtf )
            ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

            // Make reference genome
            CELLRANGER_MKREF(
                ch_fasta,
                CELLRANGER_MKGTF.out.gtf,
                "gex_reference"
            )
            ch_versions = ch_versions.mix(CELLRANGER_MKREF.out.versions)
            ch_cellranger_gex_index = CELLRANGER_MKREF.out.reference

        }
        // ch_cellranger_gex_index = PREPARE_GENOME.out.gex_index.ifEmpty { Channel.value( empty_file ) }
        // ch_cellranger_vdj_index = PREPARE_GENOME.out.vdj_index.ifEmpty { Channel.value( empty_file ) }
        // ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    emit:
        ch_versions
}
