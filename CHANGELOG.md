# nf-core/scrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- Add a checker so that `--fb_reference` does not break the pipeline in case `ab` files are not used in `cellranger multi` sub-workflow.
- Update cellbender module to latest nf-core version ([#419](https://github.com/nf-core/scrnaseq/pull/419/))
- Add profile for gpu processes ([#419](https://github.com/nf-core/scrnaseq/pull/419/))
- Update example usage command in README with valid reference genome parameter ([#420](https://github.com/nf-core/scrnaseq/issues/339))
- Document better that `cellbender` is used for empty drops calling and not the `emptydrops` method (([#339](https://github.com/nf-core/scrnaseq/issues/420)))

## v3.0.0 - 2024-12-09

## Backwards-incompatible changes

- Remove universc workflow from pipeline ([#289](https://github.com/nf-core/scrnaseq/issues/289)).
- Remove emptydrops from the pipeline, in favor of cellbender ([#369](https://github.com/nf-core/scrnaseq/pull/369)).

## Additions

- Add `--save_align_intermeds` parameter that publishes BAM files to the output directory (for `starsolo`, `cellranger` and `cellranger multi`) ([#384](https://github.com/nf-core/scrnaseq/issues/384)).

## Fixes

- Add support for pre-built indexes in `genomes.config` file for `cellranger`, `cellranger-arc`, `simpleaf` and `simpleaf txp2gene` ([#371](https://github.com/nf-core/scrnaseq/issues/371)).
- Refactor matrix conversion code. Output from all aligners is initially converted to AnnData h5ad that is used for
  downstream code such as cellbender. H5ad objects are converted to Seurat and SingleCellExperiment at the end
  using anndataR. This reduced the pipeline complexity and resolved various issues relating to output format conversion
  ([#369](https://github.com/nf-core/scrnaseq/pull/369)).
- Fix problem with `test_full` that was not running out of the box, since code was trying to overwrite parameters in the workflow, which is not possible ([#366](https://github.com/nf-core/scrnaseq/issues/366)).

## v2.7.1 - 2024-08-13

- Fix that tests have not been executed with nf-test v0.9 ([#359](https://github.com/nf-core/scrnaseq/pull/359))
- Add support for 10XV4 chemistry ([#348](https://github.com/nf-core/scrnaseq/pull/348))
- Fix issues with predefined STAR index ([#350](https://github.com/nf-core/scrnaseq/pull/350))
- Update modules ([#351](https://github.com/nf-core/scrnaseq/pull/351))
- Fix resource specifications for `cellranger mkref`/`cellrangerarc mkref` ([#352](https://github.com/nf-core/scrnaseq/pull/352))

## v2.7.0 - 2024-06-03

- Apply `check_max` to AlevinQC time limit ([#335](https://github.com/nf-core/scrnaseq/pull/335))
- Update template to v2.14.1 ([#328](https://github.com/nf-core/scrnaseq/pull/328))
- Avoid filename collisions in cellranger-arc ([#321](https://github.com/nf-core/scrnaseq/pull/321))
- Add cellranger multi subworkflow ([#247](https://github.com/nf-core/scrnaseq/issues/247))
  - Add support for 10x multiplexed, multi-omics and FFPE samples
  - Allow the use of gzipped fasta and GTF files
- Fix that pipeline couldn't run without GTF file, even when an aligner index was specified ([#322](https://github.com/nf-core/scrnaseq/pull/322))

## v2.6.0 - 2024-04-16

- Update cellranger to v8.0.0 ([#317](https://github.com/nf-core/scrnaseq/pull/317))
- Change from pytests to nf-test ([#291](https://github.com/nf-core/scrnaseq/pull/291))
- Update template to v2.13.1 ([#309](https://github.com/nf-core/scrnaseq/pull/309))
- Update to kallisto|bustools v0.28.2 ([#294](https://github.com/nf-core/scrnaseq/pull/294))
- Fix cellrangerarc matrix conversions and protocol selection ([#300](https://github.com/nf-core/scrnaseq/pull/300))
- Add new emptydrops calling module ([#301](https://github.com/nf-core/scrnaseq/pull/301))
- Update cellranger modules to latest version ([[#316](https://github.com/nf-core/scrnaseq/issues/316)])

## v2.5.1 - 2024-01-23

- Template update to v2.12 ([#298](https://github.com/nf-core/scrnaseq/pull/298)).
- Fix that cellranger workflow couldn't be run and enable CI for this workflow ([#288](https://github.com/nf-core/scrnaseq/pull/288)).
- Update modules ([#288]()https://github.com/nf-core/scrnaseq/pull/288).

## v2.5.0 - 2024-01-02

- Update template to v2.11.1 ([#279](https://github.com/nf-core/scrnaseq/pull/279))
- Add support for paired GEX+ATAC sequencing using cellranger-arc ([#274](https://github.com/nf-core/scrnaseq/pull/274))
- Increase default runtime limits for some processes ([#281](https://github.com/nf-core/scrnaseq/pull/281), [#284](https://github.com/nf-core/scrnaseq/pull/284))
- Better support for custom protocols ([#273](https://github.com/nf-core/scrnaseq/pull/273)).
  - The universc protocol is now specified via the `--protocol` flag
  - Any protocol specified is now passed to the respective aligner
  - Added a section to the documentation

## v2.4.1 - 2023-09-28

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
