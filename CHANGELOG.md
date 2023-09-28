# nf-core/scrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.4.1 - [date]

- Fix whitelist logic for dropseq ([#267](https://github.com/nf-core/scrnaseq/pull/267))
- Fix false-positive filename check in cellranger module ([#261](https://github.com/nf-core/scrnaseq/pull/261))
- Template update to v2.10 ([#269](https://github.com/nf-core/scrnaseq/pull/269))

## v2.4.0 - 2023-08-16 Lime Platinum Crab

- Fix typo causing empty version imformation for mtx_conversion subworkflow ([#254](https://github.com/nf-core/scrnaseq/pull/254))
- Add `singularity.registry = 'quay.io'` and bump NF version to 23.04.0 ([#237](https://github.com/nf-core/scrnaseq/pull/237))
- Fixed issue with file collisions while using cellranger ([#232](https://github.com/nf-core/scrnaseq/pull/232))
- Fix issue where multiqc inputs tried to access objects that did not exist ([#239](https://github.com/nf-core/scrnaseq/pull/239))
- Removed `public_aws_ecr` profile ([#242](https://github.com/nf-core/scrnaseq/pull/242))
- Include cellranger in MultiQC report ([#244](https://github.com/nf-core/scrnaseq/pull/244))
- Nf-core template update to v2.9 ([#245](https://github.com/nf-core/scrnaseq/pull/245))
- Update cellranger and fastqc module ([#246](https://github.com/nf-core/scrnaseq/pull/246)).
  The [updated cellranger module](https://github.com/nf-core/modules/pull/3537) now automatically renames input FASTQ
  files to match the expected naming conventions.

## v2.3.2 - 2023-06-07 Sepia Samarium Salmon

- Move containers for pipeline to quay.io ([#233](https://github.com/nf-core/scrnaseq/pull/233))

## v2.3.1 - 2023-06-02 Yellow Strontium Pinscher

- Add `public_aws_ecr` config for using the AWS mirror of containers where possible ([#225](https://github.com/nf-core/scrnaseq/pull/225))

## v2.3.0 Steelblue Waspaloy Dachshund

- Fix problem on samplesheet check related to amount of columns ([[#211](https://github.com/nf-core/scrnaseq/issues/211)])
- Fixed bug in starsolo output cardinality.

## v2.2.0

- Added support to output 10x count files in text format.
- Add gene symbols to count matrices
- Added UniverSC aligner to implement open-source version of Cell Ranger with wrapper for 40 technologies
- Update cellranger to v7.1.0 ([#205](https://github.com/nf-core/scrnaseq/pull/205)).

### Fixes

- Autocanceling previous CI runs when new changes are pushed.
- Fixed [#193](https://github.com/nf-core/scrnaseq/issues/177) by updating the Seurat container directive
- Fixed [#177](https://github.com/nf-core/scrnaseq/issues/177) by adjusting the channels generation and usage when skipping fastqc
- Fixed [#173](https://github.com/nf-core/scrnaseq/issues/173) by adjusting parameter type and adding them to modules.config
- Fixed [#170](https://github.com/nf-core/scrnaseq/issues/170) by adding UniverSC subworkflow using new module
- Fixed [#196](https://github.com/nf-core/scrnaseq/issues/196) by adjusting runtime requirements for AlevinQC
- Fixed [#191](https://github.com/nf-core/scrnaseq/issues/191) by updating simpleAF containers to latest version

## v2.1.0 - 2022-10-06 "Green Mercury Siberian Husky"

- Alevin workflow updated to use Alevin-Fry via simpleaf - thanks to @rob-p for supporting this and @fmalmeida implementing the support

### Fixes

- Fixed Kallistobustools workflow [#123](https://github.com/nf-core/scrnaseq/issues/123) by upgrading to nf-core/modules module
- Fixed matrix conversion error when running STAR with --soloFeatures GeneFull [#135](https://github.com/nf-core/scrnaseq/pull/135)
- Fixed seurat matrix conversion error when running with conda profile [#136](https://github.com/nf-core/scrnaseq/pull/136)
- Fixed Kallistobustools module [#116](https://github.com/nf-core/scrnaseq/issues/116). By updating nf-core module and making sure conversion modules take into account the different outputs produced by kallisto standard and non-standard workflows.
- Updated pipeline template to [nf-core/tools 2.6](https://github.com/nf-core/tools/releases/tag/2.6)

## v2.0.0 - 2022-06-17 "Gray Nickel Beagle"

- Pipeline ported to dsl2
- Template update with latest nf-core/tools v2.1
- Added cellranger v.7.0.0 subworkflow
- Added full size tests

### Fixes

- Make sure pipeline runs on multiple samples [#77](https://github.com/nf-core/scrnaseq/pull/77)
- Fix issue where STARsolo always uses 10XV2 chemistry [#60](https://github.com/nf-core/scrnaseq/issues/60)

## v1.1.0 - 2021-03-24 "Olive Mercury Corgi"

- Template update with latest nf-core/tools v1.13.2
- Parameters JSON Schema added [#42](https://github.com/nf-core/scrnaseq/issues/42)
- [25](https://github.com/nf-core/scrnaseq/issues/25) Fix small documentation error with wrong parameter for txp2gene

### Fixes

- [#20](https://github.com/nf-core/scrnaseq/issues/20) Fix Transcriptome Fasta argument not detected well
- [#21](https://github.com/nf-core/scrnaseq/issues/21) Fix `--kallisto_index` being ignored

## v1.0.0 - 2019-11-28 "Tiny Aluminium Crab"

Initial release of nf-core/scrnaseq, created with the [nf-core](http://nf-co.re/) template.
This includes the following workflow options:

- Salmon Alevin + AlevinQC
- STARSolo
- Kallisto / BUStools
