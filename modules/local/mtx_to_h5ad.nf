process MTX_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    // inputs from cellranger nf-core module does not come in a single sample dir
    // for each sample, the sub-folders and files come directly in array.
    tuple val(meta), path(inputs)
    path txp2gene
    path star_index

    output:
    path "${meta.id}/*h5ad", emit: h5ad
    path "${meta.id}/", emit: counts
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def file paths for aligners (except cellranger)
    if (params.aligner == 'kallisto') {
        mtx_matrix   = "*count/counts_unfiltered/*.mtx"
        barcodes_tsv = "*count/counts_unfiltered/*.barcodes.txt"
        features_tsv = "*count/counts_unfiltered/*.genes.txt"
    } else if (params.aligner == 'alevin') {
        mtx_matrix   = "*_alevin_results/af_quant/alevin/quants_mat.mtx"
        barcodes_tsv = "*_alevin_results/af_quant/alevin/quants_mat_rows.txt"
        features_tsv = "*_alevin_results/af_quant/alevin/quants_mat_cols.txt"
    } else if (params.aligner == 'star') {
        mtx_matrix   = "*.Solo.out/Gene*/filtered/matrix.mtx.gz"
        barcodes_tsv = "*.Solo.out/Gene*/filtered/barcodes.tsv.gz"
        features_tsv = "*.Solo.out/Gene*/filtered/features.tsv.gz"
    }

    //
    // run script
    //
    if (params.aligner == 'cellranger')
    """
    # convert file types
    mtx_to_h5ad.py \\
        --aligner ${params.aligner} \\
        --input filtered_feature_bc_matrix.h5 \\
        --sample ${meta.id} \\
        --out ${meta.id}/${meta.id}_matrix.h5ad
    """

    else if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    # convert file types
    for input_type in spliced unspliced ; do
        mtx_to_h5ad.py \\
            --aligner ${params.aligner} \\
            --sample ${meta.id} \\
            --input *count/counts_unfiltered/\${input_type}.mtx \\
            --barcode *count/counts_unfiltered/\${input_type}.barcodes.txt \\
            --feature *count/counts_unfiltered/\${input_type}.genes.txt \\
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
        --out ${meta.id}/${meta.id}_matrix.h5ad
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
