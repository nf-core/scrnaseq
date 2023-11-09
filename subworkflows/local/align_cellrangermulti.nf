// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI {
    take:
        input_data

    main:
        ch_versions = Channel.empty()

        // since we merged all data as a meta, now we have a channel per sample, which
        // every item is a meta map for each data-type
        // now we can split it back for passing as input to the module
        input_data
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

    emit:
        ch_versions
}
