# nf-core/scrnaseq: Changelog

## dev / unreleased

### Fixes

* [25](https://github.com/nf-core/scrnaseq/issues/25) Fix small documentation error with wrong parameter for txp2gene

### Dependency Updates

* Salmon to 1.1.0
* Samtools to 1.10
* GFFRead to 0.11.7
* Kallisto to 0.46.2
* BUSTools to 0.40.0

## V1.0.1 - 2019

### Fixes

* [#20](https://github.com/nf-core/scrnaseq/issues/20) Fix Transcriptome Fasta argument not detected well
* [#21](https://github.com/nf-core/scrnaseq/issues/21) Fix `--kallisto_index` being ignored

## v1.0.0 - 2019-11-28 "Tiny Aluminium Crab"

Initial release of nf-core/scrnaseq, created with the [nf-core](http://nf-co.re/) template.
This includes the following workflow options:

* Salmon Alevin + AlevinQC
* STARSolo
* Kallisto / BUStools
