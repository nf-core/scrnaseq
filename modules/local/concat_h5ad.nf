process CONCAT_H5AD {

    //
    // This module concatenates all h5ad, per type (raw, filtered, etc.) files generated during pipeline execution
    //


    tag "${meta.id}"

    label 'process_medium'

    conda "conda-forge::scanpy==1.10.2 conda-forge::python-igraph conda-forge::leidenalg"
    container "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

    input:
    tuple val(meta), path(h5ad)
    path samplesheet

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'concat_h5ad.py'

    stub:
    """
    touch combined_matrix.h5ad
    touch versions.yml
    """
}
