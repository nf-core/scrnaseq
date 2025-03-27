# nf-core/scrnaseq: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with at least 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

There is a strict requirement for the first 3 columns to match those defined in the table below.

| Column           | Description                                                                                                                                                                                                                                                                                                               |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`         | Required. Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`).                                                                                                                          |
| `fastq_1`        | Required. Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                                                                                                      |
| `fastq_2`        | Required. Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                                                                                                      |
| `expected_cells` | Optional. Number of cells expected for a sample. Must be an integer. If multiple rows are provided for the same sample, this must be the same number for all rows, i.e. the total number of expected cells for the sample.                                                                                                |
| `seq_center`     | Optional. Sequencing center for the sample. If multiple rows are provided for the same sample, this must be the same string for all rows. Samples sequenced at different centers are considered different samples and must have different identifiers. Used for STARsolo BAM outputs only. Overrides `params.seq_center`. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Expected cells

This parameter is currently supported by

- [Alevin-fry](https://alevin-fry.readthedocs.io/en/latest/generate_permit_list.html#:~:text=procedure%20described%20above.-,%2D%2Dexpect%2Dcells,-%3Cncells%3E%3A%20This)
- [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

Note that since cellranger v7, it is **not recommended** anymore to supply the `--expected-cells` parameter.

## Aligning options

By default (i.e. `--aligner simpleaf`), the pipeline uses [piscem](https://github.com/COMBINE-lab/piscem) to perform pseudo-alignment of reads to the reference genome and [Alevin-fry](https://alevin-fry.readthedocs.io/en/latest/) to perform the downstream BAM-level quantification. Then QC reports are generated with [AlevinQC](https://github.com/csoneson/alevinQC).

Other aligner options for running the pipeline are:

- [Kallisto](https://pachterlab.github.io/kallisto/about) & [Bustools](https://bustools.github.io/), where kallisto is used for alignment and bustools is used for downstream analysis
  - `--aligner kallisto`
- [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) to perform both alignment and downstream analysis.
  - `--aligner star`
- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) to perform both alignment and downstream analysis.
  - `--aligner cellranger`
- [Cellranger Multi](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi#what) to perform the alignment and downstream analysis of samples with multiple libraries at the same time using Feature Barcode technology that enables simultaneous profiling of the V(D)J repertoire, cell surface protein, antigen, and gene expression (GEX) data.

### If using cellranger

This pipeline automatically renames input FASTQ files to follow the
[naming convention by 10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input):

```
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```

For more details, see

- [this issue](https://github.com/nf-core/scrnaseq/issues/241), discussing various mechanisms to deal with non-conformant filenames
- [the README of the cellranger/count module](https://github.com/nf-core/modules/blob/master/modules/nf-core/cellranger/count/README.md) which demonstrates that renaming files does not affect the results.
- [the code for renaming files in the cellranger/count module](https://github.com/nf-core/modules/blob/master/modules/nf-core/cellranger/count/templates/cellranger_count.py)

As a sanity check, we verify that filenames of a pair of FASTQ files only differ by `R1`/`R2`.

### Support for different scRNA-seq protocols

The single-cell protocol used in the experiment can be specified using the `--protocol` flag.
For cellranger, it is recommended to stick with the default value `'auto'` for automatic detection of the protocol.
For all other aligner, you need to specify the protocol manually.

The three 10x Genomics protocols 3' v1 (`10XV1`), 3' v2 (`10XV2`), 3' v3 (`10XV3`), and 3' v4 (`10XV4`) are universally supported
by all aligners in the pipeline and mapped to the correct options automatically. If the protocol is unknown to the
nf-core pipeline, the value specified to `--protocol` is passed to the aligner _in verbatim_ to support additional protocols.

Here are some hints on running the various aligners with different protocols

#### Kallisto/bustools

The command `kb --list` shows all supported, preconfigured protocols. Additionally, a custom technology string such as
`0,0,16:0,16,26:1,0,0` can be speficied:

> Additionally kallisto bus will accept a string specifying a new technology in the format of bc:umi:seq where each of bc,umi and seq are a triplet of integers separated by a comma, denoting the file index, start and stop of the sequence used. For example to specify the 10xV2 technology we would use 0,0,16:0,16,26:1,0,0

For more details, please refer to the [Kallisto/bustools documentation](https://pachterlab.github.io/kallisto/manual#bus).

#### Simpleaf

Simpleaf has the ability to pass custom chemistries to Alevin-fry, in a slightly different format, e.g. `1{b[16]u[12]x:}2{r:}`.

For more details, see Simpleaf's paper, [He _et al._ 2023](https://doi.org/10.1093/bioinformatics/btad614) and the [detailed description](https://hackmd.io/@PI7Og0l1ReeBZu_pjQGUQQ/rJMgmvr13).

### If using cellranger-arc

#### Automatic file name detection

This pipeline currently **does not** automatically renames input FASTQ files to follow the
[naming convention by 10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input):

```
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```

Thus please make sure your files follow this naming convention.

#### Sample sheet definition

If you are using cellranger-arc you have to add the column _sample_type_ (atac for scATAC or gex for scRNA) and _fastq_barcode_ (part of the scATAC data) to your samplesheet as an input.

**Beware of the following points:**

- It is important that you give your scRNA and scATAC different [Sample Name]s.
- Check first which file is your barcode fastq file for your scATAC data ([see](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/using/fastq-input)).
- If you have more than one sequencing run then you have to give them another suffix (e.g., rep\*) to your [Sample Name] ([see](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/fastq-input#atac_quick_start)).

An example samplesheet for a dataset called test_scARC that has two sequencing runs for the scATAC and one seqeuncing run
from two lanes for the scRNA could look like this:

```csv
sample,fastq_1,fastq_2,fastq_barcode,sample_type
test_scARC,path/test_scARC_atac_rep1_S1_L001_R1_001.fastq.gz,path/test_scARC_atac_rep1_S1_L001_R2_001.fastq.gz,path/test_scARC_atac_rep1_S1_L001_I2_001.fastq.gz,atac
test_scARC,path/test_scARC_atac_rep2_S2_L001_R1_001.fastq.gz,path/test_scARC_atac_rep2_S2_L001_R2_001.fastq.gz,path/test_scARC_atac_rep2_S2_L001_I2_001.fastq.gz,atac
test_scARC,path/test_scARC_gex_S1_L001_R1_001.fastq.gz,path/test_scARC_gex_S1_L001_R2_001.fastq.gz,,gex
test_scARC,path/test_scARC_gex_S1_L002_R1_001.fastq.gz,path/test_scARC_gex_S1_L002_R2_001.fastq.gz,,gex
```

#### Config file and index

Cellranger-arc needs a reference index directory that you can provide with `--cellranger_index`.
Besure to provide the base path of the index (e.g., `--cellranger_index /PATH/TO/10X_REF/refdata-gex-GRCh38-2024-A/`).
Be aware, you can use for cellranger-arc the same index you use for cellranger ([see](https://kb.10xgenomics.com/hc/en-us/articles/4408281606797-Are-the-references-interchangeable-between-pipelines)).
Yet, a cellranger-arc index might include additional data (e.g., TF binding motifs). Therefore, please first check if
you have to create a new cellranger-arc index ([see here](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/advanced/references) for
more information)

If you decide to create a cellranger-arc index, then you need to create a config file to generate the index. The pipeline
can do this autmatically for you if you provide a `--fasta`, `--gtf`, and an optional `--motif` file. However, you can
also decide to provide your own config file with `--cellrangerarc_config`, then you also have to specify with `--cellrangerarc_reference`
the reference genome name that you have used and stated as _genome:_ in your config file.

## Running the pipeline

The minimum typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scrnaseq --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile and default `--type` and `--barcode_whitelist`. See below for more information about profiles and these options.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scrnaseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### If using Simpleaf

When building reference index, the underlying mapper of Simpleaf operates on a large number of temporary files on disk simutaneously. Therefore, for environments where CPUs and disk I/O are limited, such as on AWS, specifying `scratch=true` for the `SIMPLEAF_INDEX` module in the `conf/modules.config` file is necessary to avoid slowdowns and potential failures.

Example:

```nextflow title="conf/modules.config"
process {
  withName: SIMPLEAF_INDEX {
    scratch=true
  }
}
```

### If using cellranger-multi

#### Automatic file name detection

The pipeline is able to automatically rename input FASTQ files to follow the
[naming convention by 10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input):

```
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```

If your data already follows the expected naming convention, you can deactivate this behavior with `skip_cellranger_renaming`.

#### Sample sheet definition

If you are using cellranger-multi you have to add the column _feature_type_ to indicate which of the Feature Barcode Technology your data corresponds to:

| feature_type | description                            |
| ------------ | -------------------------------------- |
| `gex`        | Gene expression                        |
| `vdj`        | TCR/BCR profiling                      |
| `ab`         | Antibody profiling (feature barcoding) |
| `crispr`     | CRISPR capture                         |
| `cmo`        | Cell multiplexing oligos (CMO) tags    |
| `beam`       | _Currently not supported_              |

> More information on the Feature Barcode Technologies can be found here: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi

**Beware of the following points:**

- It is important that you give the same sample name for the different feature barcode technologies data that correspond to the same and should be analysed together.
- The pipeline will **automatically** generate the cellranger multi config file based on the given data.
- When working with multiplexed data (FFPE or CMO), you'll need a **second samplesheet** relating the multiplexed samples to the corresponding "physical" sample (details below). The `sample` column in the main samplesheet refers to the "physical" sample that may contain multiple multiplexed samples.

An example samplesheet could look like this:

```csv
sample,fastq_1,fastq_2,feature_type,expected_cells
PBMC_10K,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/5gex/5gex/subsampled_sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/5gex/5gex/subsampled_sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R2_001.fastq.gz,gex,1000
PBMC_10K,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/bcell/subsampled_sc5p_v2_hs_PBMC_10k_b_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/bcell/subsampled_sc5p_v2_hs_PBMC_10k_b_S1_L001_R2_001.fastq.gz,vdj,1000
PBMC_10K,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/5gex/5fb/subsampled_sc5p_v2_hs_PBMC_10k_5fb_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc/fastqs/5gex/5fb/subsampled_sc5p_v2_hs_PBMC_10k_5fb_S1_L001_R2_001.fastq.gz,ab,1000
PBMC_10K_CMO,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc_cmo/fastqs/gex_1/subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc_cmo/fastqs/gex_1/subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R2_001.fastq.gz,gex,1000
PBMC_10K_CMO,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc_cmo/fastqs/cmo/subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc_cmo/fastqs/cmo/subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz,cmo,1000
PBMC_10K_CMV,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/gex_1/subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/gex_1/subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R2_001.fastq.gz,gex,1000
PBMC_10K_CMV,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz,ab,1000
PBMC_10K_CMV,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/vdj/subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R1_001.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/vdj/subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R2_001.fastq.gz,vdj,1000
4PLEX_HUMAN,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L001_R1_001.subsampled.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L001_R2_001.subsampled.fastq.gz,gex,
4PLEX_HUMAN,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L002_R1_001.subsampled.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L002_R2_001.subsampled.fastq.gz,gex,
4PLEX_HUMAN,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L003_R1_001.subsampled.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L003_R2_001.subsampled.fastq.gz,gex,
4PLEX_HUMAN,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L004_R1_001.subsampled.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/10xgenomics/cellranger/4plex_scFFPE/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L004_R2_001.subsampled.fastq.gz,gex,
```

#### Additional samplesheet for multiplexed samples

You must provide those via a CSV with the `--cellranger_multi_barcodes` parameter. The file should look like this:

```csv
sample,multiplexed_sample_id,probe_barcode_ids,cmo_ids,description
PBMC_10K_CMO,PBMC_10K_CMO_PBMCs_human_1,,CMO301,PBMCs_human_1
PBMC_10K_CMO,PBMC_10K_CMO_PBMCs_human_2,,CMO302,PBMCs_human_2
4PLEX_HUMAN,Liver_BC1,BC001,,Healthy liver dissociated using the Miltenyi FFPE Tissue Dissociation Kit
4PLEX_HUMAN,Ovarian_BC2,BC002,,Ovarian cancer dissociated using the Miltenyi FFPE Dissociation Kit
4PLEX_HUMAN,Colorectal_BC3,BC003,,Colorectal cancer dissociated using the Miltenyi FFPE Dissociation Kit
4PLEX_HUMAN,Pancreas_BC4,BC004,,Healthy pancreas dissociated using the Miltenyi FFPE Tissue Dissociation Kit
```

The `sample` column must match the corresponding entry in the main samplesheet.

#### Additional reference data

- Cellranger multi needs a reference for **GEX and VDJ analysis**. They are calculated on the fly given the reference
  files (`--fasta` and `--gtf`) provided, but users can also provide their own with: `--cellranger_index`
  and `--cellranger_vdj_index`, for GEX and VDJ, respectively.

  > When running cellranger multi, without any VDJ data, users can also skip VDJ automated ref building with: `--skip_cellrangermulti_vdjref`.

- When working with **FFPE data**, a prob set needs to be specified via `--gex_frna_probe_set`. This file is typically
  [provided by 10x](https://www.10xgenomics.com/support/software/cell-ranger/downloads#probe-set-downloads).

- When working with **Cell Multiplexing Oligos (CMOs)**, a reference file needs to be provided via `--gex_cmo_set`. The
  default reference file, as well as a description how to write a custom one, are [available from the 10x documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi#cmo-ref)

- When working with **Feature barcoding (antibody capture)**, a reference file needs to be provided via `--fb_reference`.
  It relates each "feature" to the corresponding barcode sequence. The structure of this file is described in
  the [cellranger documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-feature-bc-analysis#feature-ref)

## Running the pipeline

The minimum typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scrnaseq --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile and default `--type` and `--barcode_whitelist`. See below for more information about profiles and these options.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scrnaseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scrnaseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scrnaseq releases page](https://github.com/nf-core/scrnaseq/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `gpu`
  - A generic configuration profile for enabling GPU execution for modules that have `process_gpu` label.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
