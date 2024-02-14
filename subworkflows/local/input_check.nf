
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

    grouped_ch =
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        // group replicate files together, modifies channel to [ val(meta), [ [reads_rep1], [reads_repN] ] ]
        .groupTuple(by: [0])

    if (params.aligner == 'cellrangerarc' ) {
        grouped_ch
        .map { meta, sample_type, sub_sample, reads -> [ meta, sample_type.flatten(), sub_sample.flatten(), reads.flatten() ] }
        .set { reads }
    } else {
        grouped_ch
        .map { meta, reads -> [ meta, reads.flatten() ] }
        .set { reads }
    }

    emit:
    reads                                     // channel: [ val(meta), [*], [ reads ] ]
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
    def sample_type      = row.sample_type ? [row.sample_type] : ['gex']

    def sub_sample = ""
    if (params.aligner == "cellrangerarc"){
        sub_sample = row.fastq_1.split("/")[-1].replaceAll("_S[0-9]+_L[0-9]+_R1_[0-9]+.fastq.gz","")
        fastqs.each{
            if(!it.name.contains(sub_sample)){
                exit 1, "ERROR: Please check input samplesheet -> Some files do not have the same sample name " +
                        "${sub_sample} in common!\n${it}"
            }
        }
    }

    fastq_meta = [ meta, fastqs ]

    if (params.aligner == "cellrangerarc"){
        fastq_meta = [ meta, sample_type, sub_sample, fastqs ]
    }

    return fastq_meta
}
