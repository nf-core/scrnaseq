process HIGHLY_VARIABLE_GENES  {
    tag "Feature selection and dimensionality reduction"
    publishDir "results/table", pattern: "*.csv", mode:'copy'
    publishDir "results/figures", pattern: "*.png", mode:'copy'

    container = 'docker.io/nfdata/sc_rnaseq:v1.0.0'

    input:
    tuple val(run_id), path(counts_filtered)

    output:
    tuple val(run_id), path("combined_matrix_DR.h5ad") , emit: h5ad
    path "UMAP_coordinates.csv", emit: UMAP
    path "UMAP_plot.png", emit: graph_UMAP
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    feature_selection_dimensionality_red.py -ad $counts_filtered

    echo "" >> versions.yml
    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
    END_VERSIONS
    feature_selection_dimensionality_red.py --version >> versions.yml

    """

    stub:
    """
    touch combined_matrix_DR.h5ad
    touch UMAP_coordinates.csv
    touch UMAP_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    feature_selection_dimensionality_red.py --version >> versions.yml

    """
}