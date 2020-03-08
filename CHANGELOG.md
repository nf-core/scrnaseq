# nf-core/scrnaseq: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

* FastQC to QC raw reads

### `Fixed`

* [25](https://github.com/nf-core/scrnaseq/issues/25) Fix small documentation error with wrong parameter for txp2gene

### `Dependencies`

* Added FastQC `0.11.9`
* Update Salmon `1.0.0` -> `1.1.0`
* Update Samtools `1.9` -> `1.10`
* Update GFFRead `0.11.6` -> `0.11.7`
* Update Kallisto `0.46.0` -> `0.46.1`
* Update BUSTools `0.39.4` -> `0.40.0`

## [1.0.1] - 2019

### `Fixed`

* [#20](https://github.com/nf-core/scrnaseq/issues/20) Fix Transcriptome Fasta argument not detected well
* [#21](https://github.com/nf-core/scrnaseq/issues/21) Fix `--kallisto_index` being ignored

## [1.0.0] - 2019-11-28 "Tiny Aluminium Crab"

Initial release of nf-core/scrnaseq, created with the [nf-core](http://nf-co.re/) template.
This includes the following workflow options:

* Salmon Alevin + AlevinQC
* STARSolo
* Kallisto / BUStools
