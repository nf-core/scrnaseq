
process AUTO_DETECT_PROTOCOL {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::jq=1.6'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6' :
        'biocontainers/jq:1.6' }"

    input:
    // the first FastQ file in `reads` is expected to contain the cell barcodes
    tuple val(meta), path(reads)
    val aligner
    path protocol_json
    path barcode_whitelist

    output:
    tuple val(meta), path(reads), emit: ch_fastq
    env PROTOCOL, emit: protocol
    env EXTRA_ARGS, emit: extra_args
    path "*.txt.gz", emit: whitelist
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    # convert protocols.json to table
    TABLE=\$(
        jq -r '
            ."$aligner" |
            to_entries[] |
            "\\(.key)\\t\\(.value.protocol//"")\\t\\(.value.whitelist//"")\\t\\(.value.extra_args//"")"
        ' "${protocol_json}"
    )

    # iterate over all protocols defined for the selected aligner
    MATCHING_FRACTIONS=\$(cut -f1 <<<"\$TABLE" | while read KEY; do
        # uncompress whitelist
        WHITELIST=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
        [ -n "\$WHITELIST" ] || continue # skip protocols without whitelist
        WHITELIST_FILE=\$(basename "\$WHITELIST")
        
        gzip -dcf "\$WHITELIST_FILE" > barcodes

        # subsample the FastQ file
        gzip -dcf "${reads[0]}" |
        awk 'FNR % 4 == 2' | # extract the read sequence from FastQ
        head -n 100000 > reads || true # the first 100k reads should suffice

        # extract the barcodes from the FastQ reads and count how many are valid barcodes
        awk -v KEY="\$KEY" -v OFS='\\t' '
            { \$0 = substr(\$0, 1, 14) } # the barcode is in the first 14 bases; 10X V2/3 barcodes are trimmed
            FILENAME == "barcodes" { barcodes[\$0] } # cache barcodes in memory
            FILENAME == "reads" && \$0 in barcodes { count++ } # count matches for each chemistry
            END { print KEY, count/FNR } # output fraction of matching barcodes for each chemitry
        ' barcodes reads

    done | sort -k2,2gr)

    # only trust the auto-detection if exactly one protocol matches
    echo -e "These were the fractions of matching barcodes by protocol:\\n\$MATCHING_FRACTIONS"
    MATCHING_PROTOCOLS_COUNT=\$(awk '\$2>=0.7' <<<"\$MATCHING_FRACTIONS" | wc -l)
    if [ \$MATCHING_PROTOCOLS_COUNT -ne 1 ]; then
         echo "ERROR: Found \$MATCHING_PROTOCOLS_COUNT matching protocols."
         exit 1
    fi
    KEY=\$(cut -f1 <<<"\$MATCHING_FRACTIONS" | head -n1)

    # extract attributes of chosen protocol
    PROTOCOL=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f2)
    WHITELIST_PATH=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
    WHITELIST=\$(basename "\$WHITELIST_PATH")
    
    # Remove all other whitelist files
    for file in \$PWD/*.txt.gz; do
        FILE_NAME=\$(basename "\$file")
        [ "\$FILE_NAME" != "\$WHITELIST" ] && rm "\$FILE_NAME"
    done

    # Copy the chosen whitelist file
    cp "\$WHITELIST" "whitelist.txt.gz"


    EXTRA_ARGS=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f4)
    echo \$PWD/*.txt.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$(jq --version | cut -d- -f2)
    END_VERSIONS
    """
}
