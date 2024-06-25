process MTX_TO_SEURAT {
    tag "$meta.id"
    label 'process_medium'

    conda "r-seurat"
    container "nf-core/seurat:4.3.0"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)

    output:
    path "${meta.id}/*.rds", emit: seuratObjects
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def aligner = params.aligner


    // Get a file to check input type. Some aligners bring arrays instead of a single file.
    def input_to_check = (inputs instanceof String) ? inputs : inputs[0]

    // check input type of inputs
    def is_emptydrops = '0'
    input_type = (input_to_check.toUriString().contains('unfiltered') || input_to_check.toUriString().contains('raw')) ? 'raw' : 'filtered'
    if ( params.aligner == 'alevin' ) { input_type = 'raw' } // alevin has its own filtering methods and mostly output a single mtx, raw here means, the base tool output
    if (input_to_check.toUriString().contains('emptydrops')) {
        input_type = 'custom_emptydrops_filter'
        is_emptydrops = '--is_emptydrops'
    }

    // def file paths for aligners. Cellranger is normally converted with the .h5 files
    // However, the emptydrops call, always generate .mtx files, thus, cellranger 'emptydrops' required a parsing
    if (params.aligner in [ "cellranger", "cellrangerarc", "cellrangermulti" ]) {

        mtx_dir  = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered/' : ''
        matrix   = "${mtx_dir}matrix.mtx*"
        barcodes = "${mtx_dir}barcodes.tsv*"
        features = "${mtx_dir}features.tsv*"

    } else if (params.aligner == 'kallisto') {

        kb_pattern = (input_type == 'raw') ? 'un' : ''
        mtx_dir    = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "counts_${kb_pattern}filtered"
        if ((input_type == 'custom_emptydrops_filter') && (params.kb_workflow != 'standard')) { mtx_dir = 'emptydrops_filtered/\${input_type}' } // dir has subdirs for non-standard workflows
        matrix     = "${mtx_dir}/*.mtx"
        barcodes   = "${mtx_dir}/*.barcodes.txt"
        features   = "${mtx_dir}/*.genes.names.txt"

        // kallisto allows the following workflows: ["standard", "lamanno", "nac"]
        // lamanno creates "spliced" and "unspliced"
        // nac creates "nascent", "ambiguous" "mature"
        // also, lamanno produces a barcodes and genes file for both spliced and unspliced
        // while nac keep only one for all the different .mtx files produced
        kb_non_standard_files = ""
        if (params.kb_workflow == "lamanno") {
            kb_non_standard_files = "spliced unspliced"
            matrix   = "${mtx_dir}/\${input_type}.mtx"
            barcodes = "${mtx_dir}/\${input_type}.barcodes.txt"
            features = "${mtx_dir}/\${input_type}.genes.txt"
        }
        if (params.kb_workflow == "nac") {
            kb_non_standard_files = "nascent ambiguous mature"
            matrix   = "${mtx_dir}/*\${input_type}.mtx"
            features = "${mtx_dir}/*.genes.txt"
        } // barcodes tsv has same pattern as standard workflow

    } else if (params.aligner == "alevin") {

        mtx_dir  = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : '*_alevin_results/af_quant/alevin'
        matrix   = "${mtx_dir}/quants_mat.mtx"
        barcodes = "${mtx_dir}/quants_mat_rows.txt"
        features = "${mtx_dir}/quants_mat_cols.txt"

    } else if (params.aligner == 'star') {

        mtx_dir  = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "${input_type}"
        suffix   = (input_type == 'custom_emptydrops_filter') ? '' : '.gz'
        matrix   = "${mtx_dir}/matrix.mtx${suffix}"
        barcodes = "${mtx_dir}/barcodes.tsv${suffix}"
        features = "${mtx_dir}/features.tsv${suffix}"

    }

    //
    // run script
    //
    """
    mkdir ${meta.id}
    """

    if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    # convert file types
    for input_type in ${kb_non_standard_files} ; do
        mtx_to_seurat.R \\
            ${matrix} \\
            ${barcodes} \\
            ${features} \\
            ${meta.id}/${meta.id}_\${input_type}_matrix.rds \\
            ${aligner} \\
            ${is_emptydrops}
    done
    """

    else
    """
    mtx_to_seurat.R \\
        $matrix \\
        $barcodes \\
        $features \\
        ${meta.id}/${meta.id}_${input_type}_matrix.rds \\
        ${aligner} \\
        ${is_emptydrops}
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.rds
    touch versions.yml
    """
}
