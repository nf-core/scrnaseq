
//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK_CELLRANGER_ARC } from '../../modules/local/samplesheet_check_cellranger_arc'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK_CELLRANGER_ARC ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_folder_channel(it) }
        .groupTuple(by: [0]) // group replicate files together, modifies channel to [ val(meta), [ [sample_ATAC], [sample_GEX] ] ]
        .map { meta, folders -> [ meta, folders.flatten() ] } // needs to flatten due to last "groupTuple", so we now have folders as a single array as expected by nf-core modules: [ val(meta), [ reads ] ]
        .set { folders }

    emit:
    folders                                                  // channel: [ val(meta), [ folders ] ]
    versions = SAMPLESHEET_CHECK_CELLRANGER_ARC.out.versions // channel: [ versions.yml ]
}


// Function to get list of [ meta, [ folder_ATAC, folder_GEX ] ]
def create_folder_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.expected_cells = row.expected_cells != null ? row.expected_cells : null
    meta.seq_center     = row.seq_center ? row.seq_center : params.seq_center

    // add path(s) of the fastq file(s) to the meta map
    def folder_meta = []
    if (!file(row.folder_ATAC).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Folder for ATAC data does not exist!\n${row.folder_ATAC}"
    }
    if (!file(row.folder_GEX).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Folder for GEX data does not exist!\n${row.folder_GEX}"
    }

    folder_meta = [ meta, [ file(row.folder_ATAC), file(row.folder_GEX) ] ]
    return folder_meta
}
