
//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK_MULTIOME {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_folder_channel(it) }
        .groupTuple(by: [0]) // group replicate files together, modifies channel to [ val(meta), [ folder_GEX, folder_ATAC ] ]
        .map { meta, folders -> [ meta, folders.flatten().unique() ] } // needs to flatten due to last "groupTuple", so we now have folder as a single array as expected by nf-core modules: [ val(meta), [ folders ] ]
        .set { folders }

    emit:
    folders                                     // channel: [ val(meta), [ folders ] ]
    versions = SAMPLESHEET_CHECK.out.versions   // channel: [ versions.yml ]
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_folder_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.single_end     = row.single_end.toBoolean()
    meta.expected_cells = row.expected_cells != null ? row.expected_cells : null
    meta.seq_center     = row.seq_center ? row.seq_center : params.seq_center

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    if (row.folder_ATAC != "" && !file(row.fastq_barcode).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Barcode FastQ (Dual index i5 read) file does not exist!\n${row.fastq_barcode}"
    }
    if (!file(row.folder_GEX).exists() && !file(row.folder_ATAC).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Either the folder for scATAC (folder_ATAC) or scRNA (folder_GEX) is missing!\n${row.folder_GEX}\n${row.folder_ATAC}"
    }

    fastq_meta = [ meta, [row.folder_GEX, row.folder_ATAC] ]
    
    return fastq_meta
}