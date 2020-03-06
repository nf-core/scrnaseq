# nf-core/scrnaseq: Usage

## Table of contents

* [nf-core/scrnaseq: Usage](#nf-corescrnaseq-usage)
  * [Table of contents](#table-of-contents)
  * [Introduction](#introduction)
  * [Running the pipeline](#running-the-pipeline)
    * [Updating the pipeline](#updating-the-pipeline)
    * [Reproducibility](#reproducibility)
  * [Main arguments](#main-arguments)
    * [`-profile`](#-profile)
    * [`--reads`](#--reads)
    * [`--aligner` (Required)](#--aligner-required)
      * [Alevin](#alevin)
        * [`--salmon_index`](#--salmon_index)
        * [`--txp2gene`](#--txp2gene)
      * [STARSolo](#starsolo)
        * [`--star_index`](#--star_index)
      * [Kallisto | BUStools](#kallisto--bustools)
        * [`--bustools_correct`](#--bustools_correct)
        * [`--skip_bustools`](#--skip_bustools)
        * [`--kallisto_gene_map`](#--kallisto_gene_map)
        * [`--kallisto_index`](#--kallisto_index)
    * [Cellular barcodes](#cellular-barcodes)
      * [`--type` to specify droplet type (Required)](#--type-to-specify-droplet-type-required)
      * [`--chemistry` (using cellranger barcodes) (Required)](#--chemistry-using-cellranger-barcodes-required)
      * [`--barcode_whitelist` for custom barcode whitelist](#--barcode_whitelist-for-custom-barcode-whitelist)
  * [Reference genomes](#reference-genomes)
    * [`--genome` (using iGenomes)](#--genome-using-igenomes)
    * [`--fasta`](#--fasta)
    * [`--gtf`](#--gtf)
    * [`--transcript_fasta`](#--transcript_fasta)
    * [`--save_reference`](#--save_reference)
    * [`--igenomes_ignore`](#--igenomes_ignore)
  * [Job resources](#job-resources)
    * [Automatic resubmission](#automatic-resubmission)
    * [Custom resource requests](#custom-resource-requests)
  * [AWS Batch specific parameters](#aws-batch-specific-parameters)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
  * [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`--email_on_fail`](#--email_on_fail)
    * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
    * [`-name`](#-name)
    * [`-resume`](#-resume)
    * [`-c`](#-c)
    * [`--custom_config_version`](#--custom_config_version)
    * [`--custom_config_base`](#--custom_config_base)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_email`](#--plaintext_email)
    * [`--monochrome_logs`](#--monochrome_logs)
    * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The minimum typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scrnaseq --reads '*_R{1,2}.fastq.gz' --fasta human.fasta --gtf human.gtf -profile docker
```

This will launch the pipeline with the `docker` configuration profile and default `--type` and `--barcode_whitelist`. See below for more information about profiles and these options.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scrnaseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scrnaseq releases page](https://github.com/nf-core/scrnaseq/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/scrnaseq`](http://hub.docker.com/r/nfcore/scrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/scrnaseq`](http://hub.docker.com/r/nfcore/scrnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--aligner` (Required)

The workflow can handle three types of methods:

* Kallisto/Bustools
* Salmon Alevin + AlevinQC
* STARsolo

To choose which one to use, please specify either `alevin`, `star` or `kallisto` as a parameter option for `--aligner`. By default, the pipeline runs the `alevin` option. Note that specifying another aligner option also requires choosing appropriate parameters (see below) for the selected option.

#### Alevin

##### `--salmon_index`

This can be used to specify a precomputed Salmon index in the pipeline, in order to skip the generation of required indices by Salmon itself.

##### `--txp2gene`

This allows the specification of a transcript to gene mapping file for Salmon Alevin and AlevinQC.

> This is not the same as the `kallisto_gene_map` parameter down below and is only used by the Salmon Alevin workflow.

#### STARSolo

##### `--star_index`

Specify a path to the precomputed STAR index.

> NB: This has to be computed with STAR Version 2.7 or later, as STARsolo was only first supported by STAR Version 2.7.

#### Kallisto | BUStools

##### `--bustools_correct`

If set to false, skip the correct steps after mapping with Kallisto.

##### `--skip_bustools`

When supplied, skip BUStools entirely.

##### `--kallisto_gene_map`

Specify a Kallisto gene mapping file here. If you don't, this will be automatically created in the Kallisto workflow when specifying a valid `--gtf` file.

##### `--kallisto_index`

Specify a path to the precomputed Kallisto index.

### Cellular barcodes

#### `--type` to specify droplet type (Required)

Currently, only 10X Genomics' chromium chemistry is supported. Drop-Seq, inDrop, etc may be supported in the future.

#### `--chemistry` (using cellranger barcodes) (Required)

To specify which chemistry (and thus barcode whitelist) to use, use the `--chemistry` flag. For example, to specify V3 chemistry (the default, as it is compatible with V2), use `--chemistry V3`.

These files were copied out of 10x Genomics' [cellranger](https://github.com/10XGenomics/cellranger) `cellranger/lib/python/cellranger/barcodes`, in some cases gzipped for simplicity across versions, and copied to `assets/whitelist`.

* V1: `737K-april-2014_rc.txt` --> gzipped --> `10x_V1_barcode_whitelist.txt.gz`
* V2: `737K-august-2016.txt` --> gzipped --> `10x_V2_barcode_whitelist.txt.gz`
* V3: `3M-february-2018.txt.gz` --> `10x_V3_barcode_whitelist.txt.gz`

#### `--barcode_whitelist` for custom barcode whitelist

If not using the 10X Genomics platform, a custom barcode whitelist can be used with `--barcode_whitelist`.

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

> Note that you need to specify either a `--genome` or `--fasta` when running the STARsolo workflow. The Kallisto and Alevin workflows can utilize a `--transcript_fasta` instead, whereas STAR needs a genomic fasta file as input in all cases.

### `--gtf`

Specify a valid GTF file for the workflow here.

### `--transcript_fasta`

If you intend to skip the generation of a transcriptomic fasta file, you can use this parameter to supply a transcriptomic fasta file here. If you don't specify this, it will be automatically generated from the supplied genomics fasta file utilizing the GTF annotation subsequently.

### `--save_reference`

Specify this parameter to save the indices created (STAR, Kallisto, Salmon) to the results.

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nfcore.slack.com/channels/scrnaseq/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
