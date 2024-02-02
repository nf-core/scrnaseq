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

    // check input type of inputs
    def is_emptydrops = '0'
    input_type = (inputs.toUriString().contains('unfiltered') || inputs.toUriString().contains('raw')) ? 'raw' : 'filtered'
    if ( params.aligner == 'alevin' ) { input_type = 'raw' } // alevin has its own filtering methods and mostly output a single mtx, raw here means, the base tool output
    if (inputs.toUriString().contains('emptydrops')) {
        input_type = 'custom_emptydrops_filter'
        is_emptydrops = '--is_emptydrops'
    }

    // def file paths for aligners. Cellranger is normally converted with the .h5 files
    // However, the emptydrops call, always generate .mtx files, thus, cellranger 'emptydrops' required a parsing
    if (params.aligner in [ 'cellranger', 'cellrangerarc' ]) {

        mtx_dir  = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "${input_type}_feature_bc_matrix"
        matrix   = "${mtx_dir}/matrix.mtx*"
        barcodes = "${mtx_dir}/barcodes.tsv*"
        features = "${mtx_dir}/features.tsv*"

    } else if (params.aligner == 'kallisto') {

        kb_pattern = (input_type == 'raw') ? 'un' : ''
        mtx_dir    = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "counts_${kb_pattern}filtered"
        matrix     = "${mtx_dir}/*.mtx"
        barcodes   = "${mtx_dir}/*.barcodes.txt"
        features   = "${mtx_dir}/*.genes.txt"

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
    for input_type in spliced unspliced ; do
        mtx_to_seurat.R \\
            *count/counts_unfiltered/\${input_type}.mtx \\
            *count/counts_unfiltered/\${input_type}.barcodes.txt \\
            *count/counts_unfiltered/\${input_type}.genes.txt \\
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
