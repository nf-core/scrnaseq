process CUSTOM_GENES {
    tag "${meta_geneset.res}_${meta_geneset.genes}"
    label 'process_medium'

    container 'docker.io/nfdata/sc_rnaseq:v1.0.1'

    input:
    tuple val(meta), path(input_h5mu)
    tuple val(meta_geneset), val(resolution), path(custom_geneset)

    output:
    path "*_features_plots.pdf", emit: feat_plot
    path "*_dotplot_r*.pdf", optional: true, emit: dotplot
    path "*_heatmap_r*.pdf", optional: true, emit: heatmap
    path "*_violin_r*.pdf", optional: true, emit: violin
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    custom_genes.py -ad ${input_h5mu} -g ${custom_geneset} -res ${resolution}

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS
    """

    stub:
    """
    touch geneset_features_plots.pdf
    touch geneset_dotplot_res.pdf
    touch geneset_heatmap_res.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustering.py --version >> versions.yml
    END_VERSIONS
    """
}
