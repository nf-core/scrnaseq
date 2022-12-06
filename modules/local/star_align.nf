process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::star=2.7.10a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10a--h9ee0642_0' :
        'quay.io/biocontainers/star:2.7.10a--h9ee0642_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path  index
    path  gtf
    path whitelist
    val protocol
    val star_feature
    val other_10x_parameters

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*.Solo.out')       , emit: counts
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_center = meta.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$meta.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def out_sam_type = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    // def read_pair = params.protocol.contains("chromium") ? "${reads[1]} ${reads[0]}" : "${reads[0]} ${reads[1]}" -- commented out to be removed is it is not being used

    // default values max percentile for UMI count 0.99 and max to min ratio for UMI count 10 taken from STARsolo usage
    def cell_filter = meta.expected_cells ? "--soloCellFilter CellRanger2.2 $meta.expected_cells 0.99 10" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reverse.join( "," )} ${forward.join( "," )} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --soloCBwhitelist <(gzip -cdf $whitelist) \\
        --soloType $protocol \\
        --soloFeatures $star_feature \\
        $other_10x_parameters \\
        $out_sam_type \\
        $ignore_gtf \\
        $seq_center \\
        $cell_filter \\
        $args \\

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    if [ -d ${prefix}.Solo.out ]; then
        # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
