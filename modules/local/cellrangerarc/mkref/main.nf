process CELLRANGERARC_MKREF {
    tag "$reference_name"
    label 'process_medium'

    container "nf-core/cellranger-arc:2.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path fasta
    path gtf
    path motifs
    path reference_config
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "config"           , emit: config
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fast_name = fasta.name
    def gtf_name = gtf.name
    def motifs_name = motifs.name
    def reference_config = reference_config.name
    def args = task.ext.args ?: ''

    if ( !reference_name ){
        reference_name = "cellrangerarc_reference"
    }

    """
    if [ $reference_config == [] ]; then
        generate_config.py \\
            --fasta $fast_name \\
            --gtf $gtf_name \\
            --motifs $motifs_name \\
            $args
    else
        if [ ! -f config ]; then
            mv -i $reference_config config
        fi
    fi

    cellranger-arc \\
        mkref \\
        --config=config \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
