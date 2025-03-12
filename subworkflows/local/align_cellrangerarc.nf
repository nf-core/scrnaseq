/*
 * Alignment with Cellranger Arc
 */

include {CELLRANGERARC_MKGTF} from "../../modules/nf-core/cellrangerarc/mkgtf/main.nf"
include {CELLRANGERARC_MKREF} from "../../modules/nf-core/cellrangerarc/mkref/main.nf"
include {CELLRANGERARC_COUNT} from "../../modules/nf-core/cellrangerarc/count/main.nf"

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
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or a bundle of a fasta file ('--fasta') + gtf file ('--gtf')."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGERARC_MKGTF( gtf )
            filtered_gtf = CELLRANGERARC_MKGTF.out.gtf
            ch_versions = ch_versions.mix(CELLRANGERARC_MKGTF.out.versions)

            // Make reference genome
            assert (( !params.cellrangerarc_reference && !cellrangerarc_config ) ||
                    ( params.cellrangerarc_reference && cellrangerarc_config ) ) :
                "If you provide a config file you also have to specific the reference name and vice versa."

            cellrangerarc_reference = 'cellrangerarc_reference'
            if ( params.cellrangerarc_reference ){
                cellrangerarc_reference = params.cellrangerarc_reference
            }

            CELLRANGERARC_MKREF( fasta, filtered_gtf, motifs, cellrangerarc_config, cellrangerarc_reference )
            ch_versions = ch_versions.mix(CELLRANGERARC_MKREF.out.versions)
            cellranger_index = CELLRANGERARC_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGERARC_COUNT (
            ch_fastq,
            cellranger_index
        )
        ch_versions = ch_versions.mix(CELLRANGERARC_COUNT.out.versions)

        // Parse the output channels to obtain filtered and raw matrices
        ch_matrices_filtered = parse_demultiplexed_output_channels( CELLRANGERARC_COUNT.out.outs, "filtered_feature_bc_matrix" )
        ch_matrices_raw      = parse_demultiplexed_output_channels( CELLRANGERARC_COUNT.out.outs, "raw_feature_bc_matrix"      )

    emit:
        ch_versions
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
        mtx_files.each{ if ( it.toString().contains("${pattern}") ) { desired_files.add( it ) } }
        [ meta_clone, desired_files ]
    }

    return out_ch
}
