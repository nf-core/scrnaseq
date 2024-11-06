//
// Include modules
//
include { CELLRANGER_MULTI as CELLRANGER_MULTI_DEMUX    } from "../../modules/nf-core/cellranger/multi/main.nf"
include { CELLRANGER_MULTI as CELLRANGER_MULTI_IMMUNE   } from "../../modules/nf-core/cellranger/multi/main.nf"
include { PARSE_CELLRANGERMULTI_SAMPLESHEET             } from "../../modules/local/parse_cellrangermulti_samplesheet.nf"
include { BAMTOFASTQ10X                                 } from '../../modules/nf-core/bamtofastq10x/main'

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_MULTI_ALIGN_VDJ {
    take:
        ch_fastq
        ch_cellranger_gex_index
        ch_cellranger_vdj_index
        ch_multi_samplesheet
        empty_file

    main:
        ch_versions    = Channel.empty()

        //
        // TODO: Include checkers for cellranger multi parameter combinations. For example, when VDJ data is given, require VDJ ref. If FFPE, require frna probe sets, etc.
        //

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

        // Add faux VDJ channel to first run cellranger without immune profiling
        ch_grouped_fastq.vdj.map { meta, fastqs ->
            def meta_clone = meta.clone()
            meta_clone.options = "[:]"
            [meta_clone, empty_file]
        }
        .first() // convert to value channel to be consumed indefinitely
        .set { ch_faux_vdj_fastq }
        // Add faux CMO channel to first run cellranger without sample demultiplexing
        ch_grouped_fastq.cmo.map { meta, fastqs ->
            def meta_clone = meta.clone()
            meta_clone.options = "[:]"
            [meta_clone, empty_file]
        }
        .first() // convert to value channel to be consumed indefinitely
        .set { ch_faux_cmo_fastq }
        // Add faux Ab channel
        ch_grouped_fastq.ab.map { meta, fastqs ->
            def meta_clone = meta.clone()
            meta_clone.options = "[:]"
            [meta_clone, empty_file]
        }
        .first() // convert to value channel to be consumed indefinitely
        .set { ch_faux_ab_fastq }

        // Assign other cellranger reference files
        ch_gex_frna_probeset      = params.gex_frna_probe_set            ? file(params.gex_frna_probe_set)            : []
        ch_gex_target_panel       = params.gex_target_panel              ? file(params.gex_target_panel)              : []
        ch_gex_cmo_set            = params.gex_cmo_set                   ? file(params.gex_cmo_set)                   : []
        ch_gex_barcodes           = params.gex_barcode_sample_assignment ? file(params.gex_barcode_sample_assignment) : []
        ch_fb_reference           = params.fb_reference                  ? file(params.fb_reference)                  : []
        ch_vdj_primer_index       = params.vdj_inner_enrichment_primers  ? file(params.vdj_inner_enrichment_primers)  : []
        ch_beam_antigen_panel_csv = [] // currently not implemented
        ch_beam_control_panel_csv = [] // currently not implemented

        // parse frna and barcode information
        if (ch_multi_samplesheet) {

            //
            // Here, we parse the received cellranger multi barcodes samplesheet.
            // We first use the get the PARSE_CELLRANGERMULTI_SAMPLESHEET module to check it and guarantee structure
            // and also split it to have one fnra/cmo .csv for each sample.
            //
            // The selection of the GEX fastqs is because samples are always expected to have at least GEX data.
            // Then, using "combined" map, which means, the "additional barcode information" of each sample, we then,
            // parse it to generate the cmo / frna samplesheets to be used by each sample.
            //
            // Here, to guarantee it and take advantage of the "FIFO"-rule and are sure that the data used in the
            // module is from the same sample from the "normal" samplesheet. We have to use the .concat().groupTuple()
            // pipe instead of .join() because .join() outputs first the arrays that could be joined and afterwards
            // the ones with "remainders", thus, we would not ensure "FIFO" and the same order.
            //
            // To guarantee this, we can define two nf-tests, one having only one sample with CMO and another with two
            // samples using CMOs, even if wrongly/repeated, but just to guarantee FIFO is working.
            //

            PARSE_CELLRANGERMULTI_SAMPLESHEET( ch_multi_samplesheet )

            ch_grouped_fastq.gex
            .map{ [it[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.cmo.flatten().map { [get_sample_id(it, "_cmo"), it ] } )
            .groupTuple()
            .map { if ( it.size() == 2 ) { it[1] } else { [] } } // a correct tuple from snippet will have: [ sample, cmo.csv ]
            .set { ch_cmo_barcode_csv }

            ch_grouped_fastq.gex
            .map{ [it[0].id] }
            .concat( PARSE_CELLRANGERMULTI_SAMPLESHEET.out.frna.flatten().map { [get_sample_id(it, "_frna"), it ] } )
            .groupTuple()
            .map { if ( it.size() == 2 ) { it[1] } else { [] } } // a correct tuple from snippet will have: [ sample, frna.csv ]
            .set { ch_frna_sample_csv }

        } else {
            ch_cmo_barcode_csv = []
            ch_frna_sample_csv = []
        }

        //
        // MODULE: cellranger multi
        //
        CELLRANGER_MULTI_DEMUX(
            ch_grouped_fastq.gex.map{ it[0] },
            ch_grouped_fastq.gex,
            ch_faux_vdj_fastq,
            ch_faux_ab_fastq,
            ch_grouped_fastq.beam,
            ch_grouped_fastq.cmo,
            ch_grouped_fastq.crispr,
            ch_cellranger_gex_index,
            ch_gex_frna_probeset,
            ch_gex_target_panel,
            ch_cellranger_vdj_index,
            ch_vdj_primer_index,
            ch_fb_reference,
            ch_beam_antigen_panel_csv,
            ch_beam_control_panel_csv,
            ch_gex_cmo_set,
            ch_cmo_barcode_csv,
            [],
            ch_frna_sample_csv,
            params.skip_cellranger_renaming
        )
        ch_versions = ch_versions.mix(CELLRANGER_MULTI_DEMUX.out.versions)
        ch_bam_files = extract_bam(CELLRANGER_MULTI_DEMUX.out.outs)

        //
        // MODULE: bam to fastq
        //
        BAMTOFASTQ10X(
            ch_bam_files
        )
        ch_versions = ch_versions.mix(BAMTOFASTQ10X.out.versions)
        ch_bamtofastq = extract_gex_fq(BAMTOFASTQ10X.out.fastq)

        ch_expanded_vdj = expand_feature_by_demultiplexed_samples(ch_grouped_fastq.vdj, ch_bamtofastq)
        ch_expanded_ab = expand_feature_by_demultiplexed_samples(ch_grouped_fastq.ab, ch_bamtofastq)
        ch_expanded_beam = expand_feature_by_demultiplexed_samples(ch_grouped_fastq.beam, ch_bamtofastq)
        ch_expanded_crispr = expand_feature_by_demultiplexed_samples(ch_grouped_fastq.crispr, ch_bamtofastq)
        
        //
        // MODULE: cellranger multi
        //
        CELLRANGER_MULTI_IMMUNE(
            ch_bamtofastq.map{ it[0] },
            ch_bamtofastq,
            ch_expanded_vdj,
            ch_expanded_ab,
            ch_expanded_beam,
            ch_faux_cmo_fastq,
            ch_expanded_crispr,
            ch_cellranger_gex_index,
            ch_gex_frna_probeset,
            ch_gex_target_panel,
            ch_cellranger_vdj_index,
            ch_vdj_primer_index,
            ch_fb_reference,
            ch_beam_antigen_panel_csv,
            ch_beam_control_panel_csv,
            [],
            [],
            [],
            [], // TODO
            params.skip_cellranger_renaming
        )
        ch_versions = ch_versions.mix(CELLRANGER_MULTI_IMMUNE.out.versions)

        //
        // Cellranger multi splits the results from each sample. So, a module execution will have: (1) a raw counts dir for all;
        // (2) a filtered counts dir PER sample; (3) a raw counts dir PER sample
        //
        // Thus, cellranger multi outputs data from all identified samples in a single channel, which will cause file collision.
        //
        // For the conversion, we should convert the resulting files of each sample, thus, now, we must parse the names
        // of the filtered 'per_sample_outs' of cellranger/multi on the split the channels raw / filtered.
        //

        // Split channels of raw and filtered to avoid file collision problems when loading the inputs in conversion modules.
        ch_matrices_filtered = parse_demultiplexed_output_channels( CELLRANGER_MULTI_IMMUNE.out.outs, "filtered_feature_bc_matrix" )
        ch_matrices_raw      = parse_demultiplexed_output_channels( CELLRANGER_MULTI_IMMUNE.out.outs, "raw_feature_bc_matrix"      )

    emit:
        ch_versions
        cellrangermulti_out = CELLRANGER_MULTI_IMMUNE.out.outs
        cellrangermulti_mtx = ch_matrices_raw.mix( ch_matrices_filtered )
}

def get_sample_id(in_ch, pattern="_cmo") {
    def bname = in_ch.baseName
    def idx = bname.lastIndexOf(pattern)
    def modified_bname = (idx != -1) ? bname[0..<idx] : bname
    return modified_bname
}

def parse_demultiplexed_output_channels(in_ch, pattern) {
    out_ch =
    in_ch.map { meta, mtx_files ->
        def desired_files = []
        mtx_files.each{ if ( it.toString().contains("${pattern}") ) { desired_files.add( it ) } }
        [ meta, desired_files ]
    }                    // separate only desired files
    .transpose()         // transpose for handling one meta/file pair at a time
    .map { meta, mtx_files ->
        def meta_clone = meta.clone()
        if ( mtx_files.toString().contains("per_sample_outs") ) {
            def demultiplexed_sample_id = mtx_files.toString().split('/per_sample_outs/')[1].split('/')[0]
            meta_clone.id = demultiplexed_sample_id.toString()
        }
        [ meta_clone, mtx_files ]
    }                    // check if output is from demultiplexed sample, if yes, correct meta.id for proper conversion naming
    .groupTuple( by: 0 ) // group it back as one file collection per sample

    return out_ch
}


def extract_bam(in_ch) {
    out_ch =
    in_ch.map { meta, bam_files ->
        def desired_files = []
        bam_files.each{ if ( it.toString().endsWith("sample_alignments.bam") ) { desired_files.add( it ) } }
        [ meta, desired_files ]
    }
    .transpose()         // transpose for handling one meta/file pair at a time
    .map { meta, bam_files ->
        def meta_clone = meta.clone()
        if ( bam_files.toString().contains("per_sample_outs") ) {
            def demux_id = bam_files.toString().split('/per_sample_outs/')[1].split('/')[0]
            meta_clone.sample_id = meta_clone.id
            meta_clone.id = demux_id.toString()
        }
        [ meta_clone, bam_files ]
    }                    // check if output is from demultiplexed sample, if yes, correct meta.id for proper conversion naming
    .groupTuple( by: 0 ) // group it back as one file collection per sample

    return out_ch
}

def extractParts(filename) {
    // convert demux dir, sample ID, lane, read, and sequence chunk to integers to sort files.
    // example: ${meta.sample_id}_0_1_XWEDGYQN/bamtofastq_S1_L002_R1_001.fastq.gz
    def matcher = filename =~ /_0_1_(\w+)\/\w+_S(\d{1})_L(\d{3})_R([12])_(\d{3})/
    if (matcher.find()) {
        def fstq_dr = matcher.group(1)
        def smpl_id = matcher.group(2).toInteger()
        def lane_id = matcher.group(3).toInteger()
        def read_id = matcher.group(4).toInteger()
        def chnk_id = matcher.group(5).toInteger()
        return [fstq_dr, smpl_id, lane_id, chnk_id, read_id]
    }
    return [0, 0, 0, 0, 0] // Default value if pattern not found
}

def extract_gex_fq(in_ch) {
    // Extract GEX fastq files from bamtofastq output and sort files in read pairs.
    def out_ch =
    in_ch.map { meta, fns ->
        def meta_clone = meta.clone()
        meta_clone.options['check-library-compatibility'] = false // in order for downstream immune profiling not to fail.
        def desired_files = []
        // GEX fq files are located in the "*_0_1_*" directory as the multi config always starts with GEX files.
        fns.each{ if ( it.toString().contains("/${meta.sample_id}_0_1_") ) { desired_files.add( it ) } }
        // Sort files to pair R1 and R2 from the same lane (L001, L002, etc) and sequence number (001, 002, etc)
        def sortedFiles = desired_files.sort { a, b ->
            def partsA = extractParts(a)
            def partsB = extractParts(b)
            return partsA[0] <=> partsB[0] ?: partsA[1] <=> partsB[1] ?: partsA[2] <=> partsB[2] ?: partsA[3] <=> partsB[3] ?: partsA[4] <=> partsB[4]
        }
        [ meta_clone, sortedFiles ]
    }

    return out_ch
}

def expand_feature_by_demultiplexed_samples(in_ch, gex_ch) {
    out_ch = 
    in_ch
        .map{ meta, fns ->
            def meta_clone = meta.clone()
            meta_clone.sample_id = meta_clone.id
            [meta_clone, fns]
        }
        .cross(gex_ch) { it[0]["sample_id"] }
        .map{ftx, gex ->
            def ftx_meta_clone = ftx[0].clone()
            def gex_meta_clone = gex[0].clone()
            assert ftx_meta_clone.sample_id == gex_meta_clone.sample_id
            ftx_meta_clone.id = gex_meta_clone.id
            [ftx_meta_clone, ftx[1]]
        }

    return out_ch
}