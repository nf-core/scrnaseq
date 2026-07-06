/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGERARC_MKGTF} from "../../../modules/nf-core/cellrangerarc/mkgtf/main.nf"
include {CELLRANGERARC_MKREF} from "../../../modules/nf-core/cellrangerarc/mkref/main.nf"
include {CELLRANGERARC_COUNT} from "../../../modules/nf-core/cellrangerarc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGERARC_ALIGN {
    take:
        fasta
        gtf
        motifs
        cellranger_index
        ch_fastq
        cellrangerarc_config

    main:
        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf')."

        if (!cellranger_index) {
            assert (( !params.cellrangerarc_reference && !cellrangerarc_config ) ||
                    ( params.cellrangerarc_reference && cellrangerarc_config ) ) :
                "If you provide a config file you also have to specific the reference name and vice versa."

            def cellrangerarc_reference = params.cellrangerarc_reference ?: 'cellrangerarc'

            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGERARC_MKGTF( gtf.map { g -> [ [ id: cellrangerarc_reference ], g ] } )
            filtered_gtf = CELLRANGERARC_MKGTF.out.gtf

            // Make reference genome (single tuple channel: meta, fasta, gtf, motifs, reference_config)
            ch_cellrangerarc_mkref = filtered_gtf
                .combine(fasta)
                .map { _mkgtf_meta, gtf_path, fasta_path ->
                    [ [ id: cellrangerarc_reference ], fasta_path, gtf_path, motifs, cellrangerarc_config ]
                }

            CELLRANGERARC_MKREF( ch_cellrangerarc_mkref )
            cellranger_index = CELLRANGERARC_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGERARC_COUNT (
            ch_fastq,
            cellranger_index
        )

        // Parse the output channels to obtain filtered and raw matrices
        ch_matrices_filtered = parse_demultiplexed_output_channels( CELLRANGERARC_COUNT.out.outs, "filtered_feature_bc_matrix" )
        ch_matrices_raw      = parse_demultiplexed_output_channels( CELLRANGERARC_COUNT.out.outs, "raw_feature_bc_matrix"      )

    emit:
        cellrangerarc_out          = CELLRANGERARC_COUNT.out.outs
        cellrangerarc_mtx_filtered = ch_matrices_filtered
        cellrangerarc_mtx_raw      = ch_matrices_raw
}

// Filter the desired files based on the pattern from an input channel
def parse_demultiplexed_output_channels(in_ch, pattern) {

    def out_ch = in_ch.map { meta, mtx_files ->
        // Set the matrix type raw/filtered in the metadata based on the pattern
        def meta_clone = meta.clone()
        meta_clone.input_type = pattern.contains('raw_') ? 'raw' : 'filtered'
        // Iterate over the matrix files and add the ones matching the pattern to the desired files list
        def desired_files = []
        mtx_files.each{ path -> if ( path.toString().contains("${pattern}") ) { desired_files.add( path ) } }
        [ meta_clone, desired_files ]
    }

    return out_ch
}
