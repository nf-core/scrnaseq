process MTX_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)
    path txp2gene
    path star_index

    output:
    tuple val(input_type), path("${meta.id}/*h5ad") , emit: h5ad
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Get a file to check input type. Some aligners bring arrays instead of a single file.
    def input_to_check = (inputs instanceof String) ? inputs : inputs[0]

    // check input type of inputs
    input_type = (input_to_check.toUriString().contains('unfiltered') || input_to_check.toUriString().contains('raw')) ? 'raw' : 'filtered'
    if ( params.aligner == 'alevin' ) { input_type = 'raw' } // alevin has its own filtering methods and mostly output a single mtx, raw here means, the base tool output
    if (input_to_check.toUriString().contains('emptydrops')) { input_type = 'custom_emptydrops_filter' }

    // def file paths for aligners. Cellranger is normally converted with the .h5 files
    // However, the emptydrops call, always generate .mtx files, thus, cellranger 'emptydrops' required a parsing
    if (params.aligner in [ 'cellranger', 'cellrangerarc', 'cellrangermulti' ] && input_type == 'custom_emptydrops_filter') {

        aligner      = 'cellranger'
        txp2gene     = ''
        star_index   = ''
        mtx_matrix   = "emptydrops_filtered/matrix.mtx"
        barcodes_tsv = "emptydrops_filtered/barcodes.tsv"
        features_tsv = "emptydrops_filtered/features.tsv"

    } else if (params.aligner == 'kallisto') {

        kb_pattern   = (input_type == 'raw') ? 'un' : ''
        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "counts_${kb_pattern}filtered"
        if ((input_type == 'custom_emptydrops_filter') && (params.kb_workflow != 'standard')) { mtx_dir = 'emptydrops_filtered/\${input_type}' } // dir has subdirs for non-standard workflows
        mtx_matrix   = "${mtx_dir}/*.mtx"
        barcodes_tsv = "${mtx_dir}/*.barcodes.txt"
        features_tsv = "${mtx_dir}/*.genes.names.txt"

        // kallisto allows the following workflows: ["standard", "lamanno", "nac"]
        // lamanno creates "spliced" and "unspliced"
        // nac creates "nascent", "ambiguous" "mature"
        // also, lamanno produces a barcodes and genes file for both spliced and unspliced
        // while nac keep only one for all the different .mtx files produced
        kb_non_standard_files = ""
        if (params.kb_workflow == "lamanno") {
            kb_non_standard_files = "spliced unspliced"
            matrix       = "${mtx_dir}/\${input_type}.mtx"
            barcodes_tsv = "${mtx_dir}/\${input_type}.barcodes.txt"
            features_tsv = "${mtx_dir}/\${input_type}.genes.txt"
        }
        if (params.kb_workflow == "nac") {
            kb_non_standard_files = "nascent ambiguous mature"
            matrix       = "${mtx_dir}/*\${input_type}.mtx"
            features_tsv = "${mtx_dir}/*.genes.txt"
        } // barcodes tsv has same pattern as standard workflow

    } else if (params.aligner == 'alevin') {

        // alevin does not have filtered/unfiltered results
        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : '*_alevin_results/af_quant/alevin'
        mtx_matrix   = "${mtx_dir}/quants_mat.mtx"
        barcodes_tsv = "${mtx_dir}/quants_mat_rows.txt"
        features_tsv = "${mtx_dir}/quants_mat_cols.txt"

    } else if (params.aligner == 'star') {

        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "${input_type}"
        suffix       = (input_type == 'custom_emptydrops_filter') ? '' : '.gz'
        mtx_matrix   = "${mtx_dir}/matrix.mtx${suffix}"
        barcodes_tsv = "${mtx_dir}/barcodes.tsv${suffix}"
        features_tsv = "${mtx_dir}/features.tsv${suffix}"

    }

    //
    // run script
    //
    if (params.aligner in [ "cellranger", "cellrangerarc", "cellrangermulti"] && input_type != 'custom_emptydrops_filter')
    """
    # convert file types
    mtx_to_h5ad.py \\
        --aligner cellranger \\
        --input *${input_type}_feature_bc_matrix.h5 \\
        --sample ${meta.id} \\
        --out ${meta.id}/${meta.id}_${input_type}_matrix.h5ad
    """

    else if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    # convert file types
    for input_type in ${kb_non_standard_files} ; do
        mtx_to_h5ad.py \\
            --aligner ${params.aligner} \\
            --sample ${meta.id} \\
            --input ${matrix} \\
            --barcode ${barcodes_tsv} \\
            --feature ${features_tsv} \\
            --txp2gene ${txp2gene} \\
            --star_index ${star_index} \\
            --out ${meta.id}/${meta.id}_\${input_type}_matrix.h5ad ;
    done
    """

    else
    """
    # convert file types
    mtx_to_h5ad.py \\
        --task_process ${task.process} \\
        --aligner ${params.aligner} \\
        --sample ${meta.id} \\
        --input $mtx_matrix \\
        --barcode $barcodes_tsv \\
        --feature $features_tsv \\
        --txp2gene ${txp2gene} \\
        --star_index ${star_index} \\
        --out ${meta.id}/${meta.id}_${input_type}_matrix.h5ad
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
