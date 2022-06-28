# nf-core/scrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.1dev

### Fixes

- Fix [#119](https://github.com/nf-core/scrnaseq/issues/119) Fix Invalid segment error using kallisto aligner

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
