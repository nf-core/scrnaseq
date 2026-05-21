/*
 * Alignment with Cellranger
 */

include { CELLRANGER_MKGTF } from "../../modules/nf-core/cellranger/mkgtf/main.nf"
include { CELLRANGER_MKREF } from "../../modules/nf-core/cellranger/mkref/main.nf"
include { CELLRANGER_COUNT } from "../../modules/nf-core/cellranger/count/main.nf"
// Modules for Velocyto launch:
include { VELOCYTO } from "../../modules/nf-core/velocyto/main.nf"
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ALIGN {
    take:
        fasta
        gtf
        cellranger_index
        ch_fastq
        protocol

    main:
        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--fasta') and a gtf file ('--gtf')."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF( gtf )

            // Make reference genome
            CELLRANGER_MKREF( fasta, CELLRANGER_MKGTF.out.gtf, "cellranger_reference" )
            cellranger_index = CELLRANGER_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGER_COUNT (
            // TODO what is `gem` and why is it needed?
            ch_fastq.map{ meta, reads -> [meta + ["chemistry": protocol, "gem": meta.id, "samples": [meta.id]], reads] },
            cellranger_index
        )

        //
        // Split channels of raw and filtered to avoid file collision problems when loading the inputs in conversion modules.
        //
        ch_matrices_raw =
        CELLRANGER_COUNT.out.outs.map { meta, mtx_files ->
            def desired_files = []
            mtx_files.each{
                if ( it.toString().contains("raw_feature_bc_matrix") ) { desired_files.add( it ) }
            }
            [ meta + [input_type: 'raw'], desired_files ]
        }

        ch_matrices_filtered =
        CELLRANGER_COUNT.out.outs.map { meta, mtx_files ->
            def desired_files = []
            mtx_files.each{
                if ( it.toString().contains("filtered_feature_bc_matrix") ) { desired_files.add( it ) }
            }
            [ meta + [input_type: 'filtered'], desired_files ]
        }

        // Run Velocyto on the output if requested by --run_velocyto true:
        if ( params.run_velocyto ) {
            // Extract the two Cell Ranger files that VELOCYTO needs:
            ch_velocyto_files =
                CELLRANGER_COUNT.out.outs
                    .map { meta, cellranger_output_files ->

                        def bam = cellranger_output_files.find {
                            it.toString().endsWith('/possorted_genome_bam.bam')
                        }

                        def barcodes = cellranger_output_files.find {
                            it.toString().endsWith('/filtered_feature_bc_matrix/barcodes.tsv.gz')
                        }

                        tuple(meta, barcodes, bam)
                    }

            // SAMTOOLS_SORT requires this extra input (index, optional).
            ch_no_index    = Channel.value('')
            ch_fasta_for_sort = fasta.map { fa -> tuple([id: 'genome'], fa) }

            // Create cellsorted_possorted_genome_bam.bam
            SAMTOOLS_SORT(
                ch_velocyto_files.map { meta, barcodes, bam -> tuple(meta, bam) },
                ch_fasta_for_sort,
                ch_no_index
            )

            // Recombine into the exact tuple VELOCYTO expects
            ch_velocyto_input =
                ch_velocyto_files
                    .join(SAMTOOLS_SORT.out.bam)
                    .map { meta, barcodes, bam, sorted_bam ->
                        tuple(meta + [input_type: 'velocyto'], barcodes, bam, sorted_bam)
                    }

            VELOCYTO(
                ch_velocyto_input,
                gtf
            )
            ch_versions = ch_versions.mix(VELOCYTO.out.versions)
        }

    emit:
        cellranger_out               = CELLRANGER_COUNT.out.outs
        cellranger_matrices_raw      = ch_matrices_raw
        cellranger_matrices_filtered = ch_matrices_filtered
        star_index                   = cellranger_index
}
