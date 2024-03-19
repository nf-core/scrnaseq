process EMPTYDROPS_CELL_CALLING {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-dropletutils"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dropletutils:1.18.0--r42hf17093f_1' :
        'quay.io/biocontainers/bioconductor-dropletutils:1.18.0--r42hf17093f_1' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("emptydrops_filtered"), emit: filtered_matrices

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.aligner == "cellranger") {

        matrix   = "matrix.mtx.gz"
        barcodes = "barcodes.tsv.gz"
        features = "features.tsv.gz"

    } else if (params.aligner == "kallisto") {

        matrix     = "counts_unfiltered/*.mtx"
        barcodes   = "counts_unfiltered/*.barcodes.txt"
        features   = "counts_unfiltered/*.genes.names.txt"

        // kallisto allows the following workflows: ["standard", "lamanno", "nac"]
        // lamanno creates "spliced" and "unspliced"
        // nac creates "nascent", "ambiguous" "mature"
        // also, lamanno produces a barcodes and genes file for both spliced and unspliced
        // while nac keep only one for all the different .mtx files produced
        kb_non_standard_files = ""
        if (params.kb_workflow == "lamanno") {
            kb_non_standard_files = "spliced unspliced"
            matrix   = "counts_unfiltered/\${input_type}.mtx"
            barcodes = "counts_unfiltered/\${input_type}.barcodes.txt"
            features = "counts_unfiltered/\${input_type}.genes.txt"
        }
        if (params.kb_workflow == "nac") {
            kb_non_standard_files = "nascent ambiguous mature"
            matrix   = "counts_unfiltered/*\${input_type}.mtx"
            features = "counts_unfiltered/*.genes.txt"
        } // barcodes tsv has same pattern as standard workflow

    } else if (params.aligner == "alevin") {

        matrix   = "*_alevin_results/af_quant/alevin/quants_mat.mtx"
        barcodes = "*_alevin_results/af_quant/alevin/quants_mat_rows.txt"
        features = "*_alevin_results/af_quant/alevin/quants_mat_cols.txt"

    } else if (params.aligner == 'star') {

        matrix   = "raw/matrix.mtx.gz"
        barcodes = "raw/barcodes.tsv.gz"
        features = "raw/features.tsv.gz"

    }

    //
    // run script
    //
    if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    # convert file types
    for input_type in ${kb_non_standard_files} ; do
        mkdir -p emptydrops_filtered/\${input_type}
        emptydrops_cell_calling.R \\
            ${matrix} \\
            ${barcodes} \\
            ${features} \\
            emptydrops_filtered/\${input_type} \\
            ${params.aligner} \\
            0
    done
    """

    else
    """
    mkdir emptydrops_filtered/
    emptydrops_cell_calling.R \\
        $matrix \\
        $barcodes \\
        $features \\
        emptydrops_filtered \\
        ${params.aligner} \\
        0
    """

    stub:
    """
    mkdir emptydrops_filtered
    touch emptydrops_filtered/empty_file
    """
}
