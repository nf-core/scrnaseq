
//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    reads = null
    versions = null

    if (params.aligner == "cellranger-arc"){
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            // group replicate files together, modifies channel to 
            // [ val(meta), [ multimeta_s1, multimeta_s1 ], [ [reads_rep1], [reads_repN] ] ]
            .groupTuple(by: [0])
            // needs to flatten due to last "groupTuple", so we now have reads as a single array as expected by 
            // nf-core modules: [ val(meta), [multi_meta], [ reads ] ]
            .map { meta, multi_meta, reads -> [ meta, multi_meta.flatten(), reads.flatten() ] } 
            .set { reads }
        versions = SAMPLESHEET_CHECK.out.versions
    } else {
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            // group replicate files together, modifies channel to [ val(meta), [ [reads_rep1], [reads_repN] ] ]
            .groupTuple(by: [0])
            // needs to flatten due to last "groupTuple", so we now have reads as a single array as expected by 
            // nf-core modules: [ val(meta), [ reads ] ]
            .map { meta, reads -> [ meta, reads.flatten() ] } 
            .set { reads }
        versions = SAMPLESHEET_CHECK.out.versions
    }

    emit:
    reads                                     // channel: [ val(meta), [multi_meta], [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}


// Function to get list of [ meta, [ multimeta ] , [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.single_end     = row.single_end.toBoolean()
    meta.expected_cells = row.expected_cells != null ? row.expected_cells : null
    meta.seq_center     = row.seq_center ? row.seq_center : params.seq_center

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    def fastqs = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastqs = [ file(row.fastq_1) ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastqs = [ file(row.fastq_1), file(row.fastq_2) ]
        if (row.sample_type == "atac") {
            if (row.fastq_barcode == "") {
                exit 1, "ERROR: Please check input samplesheet -> Barcode FastQ (Dual index i5 read) file is missing!\n"
            }
            if (!file(row.fastq_barcode).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Barcode FastQ (Dual index i5 read) file does not exist!" +
                        "\n${row.fastq_barcode}"
            }
            fastqs.add(file(row.fastq_barcode))
        }
    }

    // define meta_data for multiome
    def multi_meta  = []
    multi_meta      = row.sample_type ? [row.sample_type] : [param.sample_type]

    if (params.aligner == "cellranger-arc"){
        sub_sample = row.fastq_1.split("/")[-1].replaceAll("_S[0-9]+_L[0-9]+_R1_[0-9]+.fastq.gz","")
        fastqs.each{
            if(!it.name.contains(sub_sample)){
                exit 1, "ERROR: Please check input samplesheet -> Some files do not have the same sample name " +
                        "${sub_sample} in common!\n${it}"
            }
        }
        multi_meta.add(sub_sample)
    }

    fastq_meta = [ meta, fastqs ]

    if (params.aligner == "cellranger-arc"){
        fastq_meta = [ meta, multi_meta, fastqs ]
    }

    return fastq_meta
}