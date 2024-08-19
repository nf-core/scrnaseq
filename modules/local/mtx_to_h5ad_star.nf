process MTX_TO_H5AD_STAR {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(inputs)
    path star_index

    output:
    tuple val(meta2), path("${meta.id}/*h5ad"), emit: h5ad
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Get a file to check input type. Some aligners bring arrays instead of a single file.
    def input_to_check = (inputs instanceof String) ? inputs : inputs[0]

    // check input type of inputs
    input_type   = (input_to_check.toUriString().contains('unfiltered') || input_to_check.toUriString().contains('raw')) ? 'raw' : 'filtered'
    meta2        = meta + [input_type: input_type]

    """
    # convert file types
    mtx_to_h5ad_star.py \\
        --task_process ${task.process} \\
        --sample ${meta.id} \\
        --input ${input_type}/matrix.mtx.gz \\
        --barcode ${input_type}/barcodes.tsv.gz \\
        --feature ${input_type}/features.tsv.gz \\
        --txp2gene ${star_index}/geneInfo.tab \\
        --out ${meta.id}/${meta.id}_${input_type}_matrix.h5ad
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
