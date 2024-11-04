## v2.7.1_auto - 2024-09-23

### Fixed bugs and new features:
- Fixed that when using alevin as the aligner, if the --fasta and --gtf parameters were not explicitly provided, the pipeline would throw an error even if --genome was specified. This issue has been resolved in this update. ([366](https://github.com/nf-core/scrnaseq/issues/366))
- Fixed that MTX_TO_H5AD module cannot take ch_star_index ([370](https://github.com/nf-core/scrnaseq/issues/370))
- Add support for when filling out your pipeline parameters in Tower, Altos internal pre-built indexes for `Cell Ranger`, `Cell Ranger ARC`, `Simpleaf index`, and `Simpleaf tx2pgene` will be automatically provided if `GRCh38` or `GRCm39` genomes are provided for `--genome` parameter. There is no need to fill in the `--cellranger_index`, `--salmon_index`, and `--tx2pgene` parameters for these specific references. For any other genome references, if you want to use pre-built indexes, these parameters will still need to be filled. ([371](https://github.com/nf-core/scrnaseq/issues/371))
- Add support for two new parameters `--cat_fastq` and `save_merged_fastq` to enable sample merging.
  * Note: the pipeline can automatically merge samples without requiring these parameters to be set to true, but this **only works if the sample FastQ files have unique names**. For example, if your samplesheet defines two rows with the same id, the pipeline will usually trigger the merging process. However, if the fastq_path column contains identical FastQ file names located in different folders, the pipeline will throw an error. In such cases, these two parameters can be used to merge the FastQ files properly, and you can also save the merged files for any additional analysis.
- Add support for when we observed that under certain circumstances, samples have mixed protocols (`10XV1`, `10XV2`, and `10XV3`, `10XV4`). The original `nf-core scRNA-seq pipeline` cannot handle this situation for the `Alevin`, `Kallisto`, and `Star` aligners. Consequently, we've introduced an additional process called `AUTO_DETECT_PROTOCOL`. This process automatically detects the appropriate protocols for different samples, preventing the entire process from stalling due to mixed protocol samples being submitted.

### Known Bugs
- After completing the cellranger mkref step, not all required files are correctly copied into the cellranger_reference folder. This results in the Cell Ranger step failing to start properly due to missing essential reference files in the index. ([5564](https://support.seqera.io/support/tickets/5564))
- The `Simpleaf` version upgrade cannot be properly applied because the `Simpleaf` indexing step runs very slowly when executing the pipeline in Tower. ([5566](https://support.seqera.io/support/tickets/5566))

### Code Difference Between nf-core scrnaseq 2.7.1 and altos scrnaseq 2.7.1_auto
```
diff --git a/assets/protocols.json b/assets/protocols.json
index 0dcfdba..0552f8d 100644
--- a/assets/protocols.json
+++ b/assets/protocols.json
@@ -1,7 +1,7 @@
 {
     "alevin": {
         "10XV1": {
-            "protocol": "1{b[14]u[10]x:}2{r:}",
+            "protocol": "10xv1",
             "whitelist": "assets/whitelist/10x_V1_barcode_whitelist.txt.gz"
         },
         "10XV2": {
@@ -13,7 +13,7 @@
             "whitelist": "assets/whitelist/10x_V3_barcode_whitelist.txt.gz"
         },
         "10XV4": {
-            "protocol": "1{b[16]u[12]x:}2{r:}",
+            "protocol": "10xv4",
             "whitelist": "assets/whitelist/10x_V4_barcode_whitelist.txt.gz"
         },
         "dropseq": {
diff --git a/conf/modules.config b/conf/modules.config
index 4a9f9fe..81395a1 100644
--- a/conf/modules.config
+++ b/conf/modules.config
@@ -59,18 +59,8 @@ process {
             enabled: false
         ]
     }
-    withName: AUTO_DETECT_PROTOCOL {
-        publishDir = [
-            path: { "${params.outdir}/pipeline_info" },
-            mode: params.publish_dir_mode,
-            saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
-            enabled: false
-        ]
-    }
 }
 
-
-
 if(params.aligner == "cellranger") {
     process {
         withName: CELLRANGER_MKGTF {
diff --git a/modules/local/auto_detect_protocol.nf b/modules/local/auto_detect_protocol.nf
deleted file mode 100644
index a4cfc55..0000000
--- a/modules/local/auto_detect_protocol.nf
+++ /dev/null
@@ -1,95 +0,0 @@
-
-process AUTO_DETECT_PROTOCOL {
-    tag "$meta.id"
-    label 'process_single'
-
-    conda 'conda-forge::jq=1.6'
-    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
-        'https://depot.galaxyproject.org/singularity/jq:1.6' :
-        'biocontainers/jq:1.6' }"
-
-    input:
-    // the first FastQ file in `reads` is expected to contain the cell barcodes
-    tuple val(meta), path(reads)
-    val aligner
-    path protocol_json
-    path barcode_whitelist
-
-    output:
-    tuple val(meta), path(reads), emit: ch_fastq
-    env PROTOCOL, emit: protocol
-    env EXTRA_ARGS, emit: extra_args
-    path "*.txt.gz", emit: whitelist
-    path "versions.yml", emit: versions
-
-    when:
-    task.ext.when == null || task.ext.when
-
-    script:
-    """
-
-    # convert protocols.json to table
-    TABLE=\$(
-        jq -r '
-            ."$aligner" |
-            to_entries[] |
-            "\\(.key)\\t\\(.value.protocol//"")\\t\\(.value.whitelist//"")\\t\\(.value.extra_args//"")"
-        ' "${protocol_json}"
-    )
-
-    # iterate over all protocols defined for the selected aligner
-    MATCHING_FRACTIONS=\$(cut -f1 <<<"\$TABLE" | while read KEY; do
-        # uncompress whitelist
-        WHITELIST=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
-        [ -n "\$WHITELIST" ] || continue # skip protocols without whitelist
-        WHITELIST_FILE=\$(basename "\$WHITELIST")
-        
-        gzip -dcf "\$WHITELIST_FILE" > barcodes
-
-        # subsample the FastQ file
-        gzip -dcf "${reads[0]}" |
-        awk 'FNR % 4 == 2' | # extract the read sequence from FastQ
-        head -n 100000 > reads || true # the first 100k reads should suffice
-
-        # extract the barcodes from the FastQ reads and count how many are valid barcodes
-        awk -v KEY="\$KEY" -v OFS='\\t' '
-            { \$0 = substr(\$0, 1, 14) } # the barcode is in the first 14 bases; 10X V2/3 barcodes are trimmed
-            FILENAME == "barcodes" { barcodes[\$0] } # cache barcodes in memory
-            FILENAME == "reads" && \$0 in barcodes { count++ } # count matches for each chemistry
-            END { print KEY, count/FNR } # output fraction of matching barcodes for each chemitry
-        ' barcodes reads
-
-    done | sort -k2,2gr)
-
-    # only trust the auto-detection if exactly one protocol matches
-    echo -e "These were the fractions of matching barcodes by protocol:\\n\$MATCHING_FRACTIONS"
-    MATCHING_PROTOCOLS_COUNT=\$(awk '\$2>=0.7' <<<"\$MATCHING_FRACTIONS" | wc -l)
-    if [ \$MATCHING_PROTOCOLS_COUNT -ne 1 ]; then
-         echo "ERROR: Found \$MATCHING_PROTOCOLS_COUNT matching protocols."
-         exit 1
-    fi
-    KEY=\$(cut -f1 <<<"\$MATCHING_FRACTIONS" | head -n1)
-
-    # extract attributes of chosen protocol
-    PROTOCOL=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f2)
-    WHITELIST_PATH=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
-    WHITELIST=\$(basename "\$WHITELIST_PATH")
-    
-    # Remove all other whitelist files
-    for file in \$PWD/*.txt.gz; do
-        FILE_NAME=\$(basename "\$file")
-        [ "\$FILE_NAME" != "\$WHITELIST" ] && rm "\$FILE_NAME"
-    done
-
-    # Copy the chosen whitelist file
-    cp "\$WHITELIST" "whitelist.txt.gz"
-
-
-    EXTRA_ARGS=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f4)
-    echo \$PWD/*.txt.gz
-    cat <<-END_VERSIONS > versions.yml
-    "${task.process}":
-        jq: \$(jq --version | cut -d- -f2)
-    END_VERSIONS
-    """
-}
diff --git a/modules/local/simpleaf_index.nf b/modules/local/simpleaf_index.nf
index 139b726..8e8bd51 100644
--- a/modules/local/simpleaf_index.nf
+++ b/modules/local/simpleaf_index.nf
@@ -7,7 +7,6 @@ process SIMPLEAF_INDEX {
         'https://depot.galaxyproject.org/singularity/simpleaf:0.10.0--h9f5acd7_1' :
         'biocontainers/simpleaf:0.10.0--h9f5acd7_1' }"
 
-
     input:
     path genome_fasta
     path transcript_fasta
@@ -34,7 +33,6 @@ process SIMPLEAF_INDEX {
     simpleaf set-paths
 
     # run simpleaf index
-
     simpleaf \\
         index \\
         --threads $task.cpus \\
diff --git a/modules/nf-core/cat/fastq/environment.yml b/modules/nf-core/cat/fastq/environment.yml
deleted file mode 100644
index 8c69b12..0000000
--- a/modules/nf-core/cat/fastq/environment.yml
+++ /dev/null
@@ -1,7 +0,0 @@
-name: cat_fastq
-channels:
-  - conda-forge
-  - bioconda
-  - defaults
-dependencies:
-  - conda-forge::coreutils=8.30
diff --git a/modules/nf-core/cat/fastq/main.nf b/modules/nf-core/cat/fastq/main.nf
deleted file mode 100644
index 4b3dc1b..0000000
--- a/modules/nf-core/cat/fastq/main.nf
+++ /dev/null
@@ -1,84 +0,0 @@
-process CAT_FASTQ {
-    tag "$meta.id"
-    label 'process_single'
-
-    conda "${moduleDir}/environment.yml"
-    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
-        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
-        'nf-core/ubuntu:20.04' }"
-
-    input:
-    tuple val(meta), path(reads, stageAs: "input*/*")
-
-    output:
-    tuple val(meta), path("*_merged_S1_R*.fastq.gz"), emit: reads
-    path "versions.yml"                       , emit: versions
-
-    when:
-    task.ext.when == null || task.ext.when
-
-    script:
-    def args = task.ext.args ?: ''
-    def prefix = task.ext.prefix ?: "${meta.id}"
-    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
-    if (meta.single_end) {
-        if (readList.size >= 1) {
-            """
-            cat ${readList.join(' ')} > ${prefix}_merged_S1_R1_001.fastq.gz
-
-            cat <<-END_VERSIONS > versions.yml
-            "${task.process}":
-                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
-            END_VERSIONS
-            """
-        }
-    } else {
-        if (readList.size >= 2) {
-            def read1 = []
-            def read2 = []
-            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
-            """
-            cat ${read1.join(' ')} > ${prefix}_merged_S1_R1_001.fastq.gz
-            cat ${read2.join(' ')} > ${prefix}_merged_S1_R2_001.fastq.gz
-
-            cat <<-END_VERSIONS > versions.yml
-            "${task.process}":
-                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
-            END_VERSIONS
-            """
-        }
-    }
-
-    stub:
-    def prefix = task.ext.prefix ?: "${meta.id}"
-    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
-    
-    """
-    echo "Prefix: ${prefix}"
-    """
-    
-    if (meta.single_end) {
-        if (readList.size >= 1) {
-            """
-            echo '' | gzip > ${prefix}_merged_S1_R1_001.fastq.gz
-
-            cat <<-END_VERSIONS > versions.yml
-            "${task.process}":
-                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
-            END_VERSIONS
-            """
-        }
-    } else {
-        if (readList.size >= 2) {
-            """
-            echo '' | gzip > ${prefix}_merged_S1_R1_001.fastq.gz
-            echo '' | gzip > ${prefix}_merged_S1_R2_001.fastq.gz
-
-            cat <<-END_VERSIONS > versions.yml
-            "${task.process}":
-                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
-            END_VERSIONS
-            """
-        }
-    }
-}
diff --git a/modules/nf-core/cat/fastq/meta.yml b/modules/nf-core/cat/fastq/meta.yml
deleted file mode 100644
index db4ac3c..0000000
--- a/modules/nf-core/cat/fastq/meta.yml
+++ /dev/null
@@ -1,42 +0,0 @@
-name: cat_fastq
-description: Concatenates fastq files
-keywords:
-  - cat
-  - fastq
-  - concatenate
-tools:
-  - cat:
-      description: |
-        The cat utility reads files sequentially, writing them to the standard output.
-      documentation: https://www.gnu.org/software/coreutils/manual/html_node/cat-invocation.html
-      licence: ["GPL-3.0-or-later"]
-input:
-  - meta:
-      type: map
-      description: |
-        Groovy Map containing sample information
-        e.g. [ id:'test', single_end:false ]
-  - reads:
-      type: file
-      description: |
-        List of input FastQ files to be concatenated.
-output:
-  - meta:
-      type: map
-      description: |
-        Groovy Map containing sample information
-        e.g. [ id:'test', single_end:false ]
-  - reads:
-      type: file
-      description: Merged fastq file
-      pattern: "*.{merged.fastq.gz}"
-  - versions:
-      type: file
-      description: File containing software versions
-      pattern: "versions.yml"
-authors:
-  - "@joseespinosa"
-  - "@drpatelh"
-maintainers:
-  - "@joseespinosa"
-  - "@drpatelh"
diff --git a/modules/nf-core/cat/fastq/nextflow.config b/modules/nf-core/cat/fastq/nextflow.config
deleted file mode 100644
index 64de2e4..0000000
--- a/modules/nf-core/cat/fastq/nextflow.config
+++ /dev/null
@@ -1,9 +0,0 @@
-process {
-    withName: 'CAT_FASTQ' {
-        publishDir = [
-            path: { params.save_merged_fastq ? "${params.outdir}/fastq" : params.outdir },
-            mode: params.publish_dir_mode,
-            saveAs: { filename -> (filename.endsWith('.fastq.gz') && params.save_merged_fastq) ? filename : null }
-        ]
-    }
-}
diff --git a/modules/nf-core/cat/fastq/tests/main.nf.test b/modules/nf-core/cat/fastq/tests/main.nf.test
deleted file mode 100644
index 6cc7aad..0000000
--- a/modules/nf-core/cat/fastq/tests/main.nf.test
+++ /dev/null
@@ -1,244 +0,0 @@
-// NOTE The version snaps may not be consistant
-// https://github.com/nf-core/modules/pull/4087#issuecomment-1767948035
-nextflow_process {
-
-    name "Test Process CAT_FASTQ"
-    script "../main.nf"
-    process "CAT_FASTQ"
-
-    test("test_cat_fastq_single_end") {
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_paired_end") {
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:false ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_single_end_same_name") {
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_paired_end_same_name") {
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:false ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_single_end_single_file") {
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_single_end - stub") {
-
-        options "-stub"
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_paired_end - stub") {
-
-        options "-stub"
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:false ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_single_end_same_name - stub") {
-
-        options "-stub"
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_paired_end_same_name - stub") {
-
-        options "-stub"
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:false ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
-                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-
-    test("test_cat_fastq_single_end_single_file - stub") {
-
-        options "-stub"
-
-        when {
-            process {
-                """
-                input[0] = Channel.of([
-                    [ id:'test', single_end:true ], // meta map
-                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)]
-                ])
-                """
-            }
-        }
-
-        then {
-            assertAll(
-                { assert process.success },
-                { assert snapshot(process.out).match() }
-            )
-        }
-    }
-}
diff --git a/modules/nf-core/cat/fastq/tests/main.nf.test.snap b/modules/nf-core/cat/fastq/tests/main.nf.test.snap
deleted file mode 100644
index aec119a..0000000
--- a/modules/nf-core/cat/fastq/tests/main.nf.test.snap
+++ /dev/null
@@ -1,376 +0,0 @@
-{
-    "test_cat_fastq_single_end": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,ee314a9bd568d06617171b0c85f508da"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,ee314a9bd568d06617171b0c85f508da"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-01-17T17:30:39.816981"
-    },
-    "test_cat_fastq_single_end_same_name": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-01-17T17:32:35.229332"
-    },
-    "test_cat_fastq_single_end_single_file": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,4161df271f9bfcd25d5845a1e220dbec"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,4161df271f9bfcd25d5845a1e220dbec"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-01-17T17:34:00.058829"
-    },
-    "test_cat_fastq_paired_end_same_name": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22",
-                            "test_2.merged.fastq.gz:md5,a52cab0b840c7178b0ea83df1fdbe8d5"
-                        ]
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22",
-                            "test_2.merged.fastq.gz:md5,a52cab0b840c7178b0ea83df1fdbe8d5"
-                        ]
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-01-17T17:33:33.031555"
-    },
-    "test_cat_fastq_single_end - stub": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-07-05T12:07:28.244999"
-    },
-    "test_cat_fastq_paired_end_same_name - stub": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940",
-                            "test_2.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                        ]
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940",
-                            "test_2.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                        ]
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-07-05T12:07:57.070911"
-    },
-    "test_cat_fastq_single_end_same_name - stub": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-07-05T12:07:46.796254"
-    },
-    "test_cat_fastq_paired_end": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22",
-                            "test_2.merged.fastq.gz:md5,a52cab0b840c7178b0ea83df1fdbe8d5"
-                        ]
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,3ad9406595fafec8172368f9cd0b6a22",
-                            "test_2.merged.fastq.gz:md5,a52cab0b840c7178b0ea83df1fdbe8d5"
-                        ]
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-01-17T17:32:02.270935"
-    },
-    "test_cat_fastq_paired_end - stub": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940",
-                            "test_2.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                        ]
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": false
-                        },
-                        [
-                            "test_1.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940",
-                            "test_2.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                        ]
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-07-05T12:07:37.807553"
-    },
-    "test_cat_fastq_single_end_single_file - stub": {
-        "content": [
-            {
-                "0": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "1": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ],
-                "reads": [
-                    [
-                        {
-                            "id": "test",
-                            "single_end": true
-                        },
-                        "test.merged.fastq.gz:md5,68b329da9893e34099c7d8ad5cb9c940"
-                    ]
-                ],
-                "versions": [
-                    "versions.yml:md5,d42d6e24d67004608495883e00bd501b"
-                ]
-            }
-        ],
-        "meta": {
-            "nf-test": "0.8.4",
-            "nextflow": "24.04.2"
-        },
-        "timestamp": "2024-07-05T12:14:51.861264"
-    }
-}
\ No newline at end of file
diff --git a/modules/nf-core/cellranger/mkref/main.nf b/modules/nf-core/cellranger/mkref/main.nf
index d793170..a719b77 100644
--- a/modules/nf-core/cellranger/mkref/main.nf
+++ b/modules/nf-core/cellranger/mkref/main.nf
@@ -35,7 +35,6 @@ process CELLRANGER_MKREF {
         --localmem=${task.memory.toGiga()} \\
         --nthreads=${task.cpus} \\
         $args
-    
 
     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
diff --git a/nextflow.config b/nextflow.config
index 7f9491d..e1b608d 100644
--- a/nextflow.config
+++ b/nextflow.config
@@ -14,20 +14,18 @@ params {
     input             = null
     save_reference    = false
     protocol          = 'auto'
-    cat_fastq         = false
-    save_merged_fastq = false
 
     // reference files
     genome            = null
     transcript_fasta  = null
-    // txp2gene          = null
-    // fasta             = null
-    // gtf               = null
+    txp2gene          = null
+    fasta             = null
+    gtf               = null
 
     // salmon alevin parameters (simpleaf)
     simpleaf_rlen     = 91
     barcode_whitelist = null
-    // salmon_index      = null
+    salmon_index      = null
 
     // kallisto bustools parameters
     kallisto_index    = null
@@ -37,14 +35,14 @@ params {
     kb_filter         = false
 
     // STARsolo parameters
-    // star_index          = null
+    star_index          = null
     star_ignore_sjdbgtf = null
     seq_center          = null
     star_feature        = "Gene"
 
     // Cellranger parameters
     skip_cellranger_renaming = false
-    // cellranger_index         = null
+    cellranger_index         = null
 
     // Cellranger ARC parameters
     motifs                  = null
diff --git a/nextflow_schema.json b/nextflow_schema.json
index 79864a4..e5fb71b 100644
--- a/nextflow_schema.json
+++ b/nextflow_schema.json
@@ -28,18 +28,6 @@
                     "description": "The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
                     "fa_icon": "fas fa-folder-open"
                 },
-                "cat_fastq": {
-                    "type": "boolean",
-                    "fa_icon": "fas fa-cat",
-                    "default": false,
-                    "description": "Concatenate FastQ files before processing."
-
-                },
-                "save_merged_fastq": {
-                    "type": "boolean",
-                    "fa_icon": "fas fa-save",
-                    "description": "Save FastQ files after merging re-sequenced libraries in the results directory."
-                },
                 "email": {
                     "type": "string",
                     "description": "Email address for completion summary.",
@@ -77,8 +65,8 @@
                 },
                 "protocol": {
                     "type": "string",
-                    "description": "The protocol used to generate the single cell data, e.g., 10x Genomics v2 Chemistry.\n\nThe protocol can be set to 'auto' (for cellranger, cellranger arc, and cellranger multi) or manually specified as '10XV1', '10XV2', '10XV3', '10XV4', or any other protocol string which will be directly passed to the respective aligner.\n\nWith the introduction of the Chemistry Auto Detection module, the pipeline now automatically detects the chemistry when using the aligners Alevin, STAR, and Kallisto. For these aligners, you can leave the protocol as 'auto', and the module will determine the appropriate chemistry automatically. However, for other cases such as cellranger, cellranger arc, and cellranger multi, the Chemistry Auto Detection module will not be triggered, and the protocol you specify (including 'auto') will behave as per the default behavior of those aligners (i.e., they will calculate the chemistry internally if set to 'auto').\n\nFor other aligners, any protocol value that is not recognized will be passed verbatim to the aligner, supporting different sequencing platforms. For more details, refer to the documentation for [kallisto](https://pachterlab.github.io/kallisto/manual#bus), [simpleaf](https://simpleaf.readthedocs.io/en/latest/quant-command.html#a-note-on-the-chemistry-flag), [starsolo](https://gensoft.pasteur.fr/docs/STAR/2.7.9a/STARsolo.html), and [universc](https://github.com/minoda-lab/universc#pre-set-configurations).",
-                    "help_text": "The default behavior is to auto-detect the protocol when running aligners like Alevin, STAR, and Kallisto, using the Chemistry Auto Detection module. For these aligners, the protocol does not need to be specified manually when set to 'auto'.\n\nFor other aligners, such as cellranger, cellranger arc, and cellranger multi, the Chemistry Auto Detection module will not be activated, and the aligner's own internal detection will apply if 'auto' is used. Alternatively, the protocol can be manually specified, and the following protocols are recognized and mapped: '10XV1', '10XV2', '10XV3', '10XV4'.\n\nAny other protocol value is passed directly to the aligner to accommodate other sequencing platforms.",
+                    "description": "The protocol that was used to generate the single cell data, e.g. 10x Genomics v2 Chemistry.\n\n Can be 'auto' (cellranger only), '10XV1', '10XV2', '10XV3', '10XV4', or any other protocol string that will get directly passed the respective aligner.",
+                    "help_text": "The default is to auto-detect the protocol when running cellranger. For all other aligners the protocol MUST be manually specified. \n\n The following protocols are recognized by the pipeline and mapped to the corresponding protocol name of the respective aligner: '10XV1', '10XV2', '10XV3', '10XV4'. \n\nAny other protocol value is passed to the aligner in verbatim to support other sequencing platforms. See the [kallisto](https://pachterlab.github.io/kallisto/manual#bus), [simpleaf](https://simpleaf.readthedocs.io/en/latest/quant-command.html#a-note-on-the-chemistry-flag), [starsolo](https://gensoft.pasteur.fr/docs/STAR/2.7.9a/STARsolo.html), and [universc](https://github.com/minoda-lab/universc#pre-set-configurations) documentations for more details.",
                     "default": "auto",
                     "fa_icon": "fas fa-cogs"
                 }
diff --git a/subworkflows/local/utils_nfcore_scrnaseq_pipeline/main.nf b/subworkflows/local/utils_nfcore_scrnaseq_pipeline/main.nf
index 1e9ac1e..dd54d35 100644
--- a/subworkflows/local/utils_nfcore_scrnaseq_pipeline/main.nf
+++ b/subworkflows/local/utils_nfcore_scrnaseq_pipeline/main.nf
@@ -207,9 +207,6 @@ def getGenomeAttribute(attribute) {
     return null
 }
 
-
-
-
 //
 // Exit pipeline if incorrect --genome key provided
 //
diff --git a/workflows/scrnaseq.nf b/workflows/scrnaseq.nf
index b511ce5..10ced22 100644
--- a/workflows/scrnaseq.nf
+++ b/workflows/scrnaseq.nf
@@ -17,45 +17,25 @@ include { softwareVersionsToYAML             } from '../subworkflows/nf-core/uti
 include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
 include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
 include { getGenomeAttribute                 } from '../subworkflows/local/utils_nfcore_scrnaseq_pipeline'
-include { AUTO_DETECT_PROTOCOL               } from '../modules/local/auto_detect_protocol'
-include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main'
+
 
 
 workflow SCRNASEQ {
 
     take:
-    ch_reads
+    ch_fastq
 
     main:
-    
-    if (params.cat_fastq){
-        ch_reads
-        .branch {
-            meta, fastqs ->
-                single  : fastqs.size() == 2
-                    return [ meta, fastqs.flatten() ]
-                multiple: fastqs.size() > 2
-                    return [ meta, fastqs.flatten() ]
-        }
-        .set { ch_fastq }
-    } else {
-        ch_fastq = ch_reads
-    }
-    
 
     protocol_config = Utils.getProtocol(workflow, log, params.aligner, params.protocol)
+    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
+        error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
+    }
 
-
-   
     params.fasta = getGenomeAttribute('fasta')
     params.gtf = getGenomeAttribute('gtf')
     params.star_index = getGenomeAttribute('star')
-    
-    // Add alevin index and txp2gene file
-    params.salmon_index = getGenomeAttribute('simpleaf')
-    params.txp2gene = getGenomeAttribute('simpleaf_tx2pgene')
 
-    
     ch_genome_fasta = params.fasta ? file(params.fasta, checkIfExists: true) : []
     ch_gtf = params.gtf ? file(params.gtf, checkIfExists: true) : []
 
@@ -93,15 +73,9 @@ workflow SCRNASEQ {
     //star params
     star_index = params.star_index ? file(params.star_index, checkIfExists: true) : null
     ch_star_index = star_index ? [[id: star_index.baseName], star_index] : []
-
     star_feature = params.star_feature
 
     //cellranger params
-    if (params.aligner in ["cellranger", "cellrangermulti"]){
-        params.cellranger_index = getGenomeAttribute('cellranger')
-    } else if (params.aligner == "cellrangerarc"){
-        params.cellranger_index = getGenomeAttribute('cellrangerarc')
-    }
     ch_cellranger_index = params.cellranger_index ? file(params.cellranger_index) : []
 
     //universc params
@@ -115,20 +89,6 @@ workflow SCRNASEQ {
     ch_versions     = Channel.empty()
     ch_mtx_matrices = Channel.empty()
 
-    if (params.cat_fastq) {
-        CAT_FASTQ (
-            ch_fastq.multiple
-        )
-        .reads
-        .mix(ch_fastq.single)
-        .set { ch_fastq }
-    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
-
-    }
-    
-
-    
-
     // Run FastQC
     if (!params.skip_fastqc) {
         FASTQC_CHECK ( ch_fastq )
@@ -136,21 +96,6 @@ workflow SCRNASEQ {
         ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CHECK.out.fastqc_multiqc.flatten())
     }
 
-
-    if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti", "universc"]) {
-        ch_protocol_json   = file( "$projectDir/assets/protocols.json", checkIfExists: true )
-        ch_whitlist = file ("$projectDir/assets/whitelist/*.gz", checkIfExists: true)
-
-        AUTO_DETECT_PROTOCOL(ch_fastq, params.aligner, ch_protocol_json, ch_whitlist)
-        protocol_config['protocol'] = AUTO_DETECT_PROTOCOL.out.protocol
-        ch_fastq = AUTO_DETECT_PROTOCOL.out.ch_fastq
-        ch_barcode_whitelist = AUTO_DETECT_PROTOCOL.out.whitelist
-        ch_extra_args = AUTO_DETECT_PROTOCOL.out.extra_args
-       
-    } else if (protocol_config['protocol'] == 'auto' && params.aligner == "universc") {
-        error "Auto-detection of protocol is not supported for UniversC"
-    }
-
     //
     // Uncompress genome fasta file if required
     //
@@ -371,7 +316,7 @@ workflow SCRNASEQ {
         ch_mtx_matrices,
         ch_input,
         ch_txp2gene,
-        star_index
+        ch_star_index
     )
 
     //Add Versions from MTX Conversion workflow too
```