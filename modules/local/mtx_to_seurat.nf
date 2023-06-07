process MTX_TO_SEURAT {
    tag "$meta.id"
    label 'process_medium'

    conda "r-seurat"
    container 'docker.io/satijalab/seurat:4.3.0'

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
    if (params.aligner == "cellranger") {
        matrix   = "matrix.mtx.gz"
        barcodes = "barcodes.tsv.gz"
        features = "features.tsv.gz"
    } else if (params.aligner == "kallisto") {
        matrix   = "*count/counts_unfiltered/*.mtx"
        barcodes = "*count/counts_unfiltered/*.barcodes.txt"
        features = "*count/counts_unfiltered/*.genes.txt"
    } else if (params.aligner == "alevin") {
        matrix   = "*_alevin_results/af_quant/alevin/quants_mat.mtx"
        barcodes = "*_alevin_results/af_quant/alevin/quants_mat_rows.txt"
        features = "*_alevin_results/af_quant/alevin/quants_mat_cols.txt"
    } else if (params.aligner == 'star') {
        matrix   = "*.Solo.out/Gene*/filtered/matrix.mtx.gz"
        barcodes = "*.Solo.out/Gene*/filtered/barcodes.tsv.gz"
        features = "*.Solo.out/Gene*/filtered/features.tsv.gz"
    }
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
            ${aligner}
    done
    """

    else
    """
    mtx_to_seurat.R \\
        $matrix \\
        $barcodes \\
        $features \\
        ${meta.id}/${meta.id}_matrix.rds \\
        ${aligner}
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.rds
    touch versions.yml
    """
}
