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

        matrix   = "raw_feature_bc_matrix/matrix.mtx.gz"
        barcodes = "raw_feature_bc_matrix/barcodes.tsv.gz"
        features = "raw_feature_bc_matrix/features.tsv.gz"

    } else if (params.aligner == "kallisto") {

        matrix     = "counts_unfiltered/*.mtx"
        barcodes   = "counts_unfiltered/*.barcodes.txt"
        features   = "counts_unfiltered/*.genes.txt"

    } else if (params.aligner == "alevin") {

        matrix   = "*_alevin_results/af_quant/alevin/quants_mat.mtx"
        barcodes = "*_alevin_results/af_quant/alevin/quants_mat_rows.txt"
        features = "*_alevin_results/af_quant/alevin/quants_mat_cols.txt"

    } else if (params.aligner == 'star') {

        matrix   = "*.Solo.out/Gene*/raw/matrix.mtx.gz"
        barcodes = "*.Solo.out/Gene*/raw/barcodes.tsv.gz"
        features = "*.Solo.out/Gene*/raw/features.tsv.gz"

    }

    //
    // run script
    //
    if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    mkdir emptydrops_filtered/
    # convert file types
    for splice_type in spliced unspliced ; do
        emptydrops_cell_calling.R \\
            *count/counts_unfiltered/\${splice_type}.mtx \\
            *count/counts_unfiltered/\${splice_type}.barcodes.txt \\
            *count/counts_unfiltered/\${splice_type}.genes.txt \\
            emptydrops_filtered \\
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
    touch emptydrops_filtered/*
    """
}
