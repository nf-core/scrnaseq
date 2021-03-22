# ![nf-core/scrnaseq](docs/images/nf-core-scrnaseq_logo.png)

**A fully automated Nextflow pipeline for Droplet-based (e.g. 10x Genomics) single-cell RNA-Seq data**.

[![GitHub Actions CI Status](https://github.com/nf-core/scrnaseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/scrnaseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/scrnaseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/scrnaseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/scrnaseq.svg)](https://hub.docker.com/r/nfcore/scrnaseq)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23scrnaseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/scrnaseq)

[![Join us on Slack](https://img.shields.io/badge/slack-nfcore/scrnaseq-blue.svg)](https://nfcore.slack.com/channels/scrnaseq)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/scrnaseq** is a bioinformatics best-practise analysis pipeline for

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

This is a community effort in building a pipeline capable to support:

* Alevin + AlevinQC
* STARSolo
* Kallisto + BUStools

## Pipeline Summary

By default, the pipeline currently performs the following:

<!-- TODO nf-core: Fill in short bullet-pointed list of default steps of pipeline -->

* Sequencing quality control (`FastQC`)
* Overall pipeline run summaries (`MultiQC`)

## Documentation

The nf-core/scrnaseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/scrnaseq/usage) and [output](https://nf-co.re/scrnaseq/output).

## Credits

The `nf-core/scrnaseq` was initiated by [Peter J. Bailey](https://github.com/PeterBailey) (Salmon Alevin, AlevinQC) with major contributions from [Olga Botvinnik](https://github.com/olgabot) (STARsolo, Testdata) and [Alex Peltzer](https://github.com/apeltzer) (Kallisto/BusTools workflow).

We thank the following people for their extensive assistance in the development
of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#scrnaseq` channel](https://nfcore.slack.com/channels/scrnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

The basic benchmarks that were used as motivation for incorporating the three available modular workflows can be found in [this publication](https://www.biorxiv.org/content/10.1101/673285v2).

We offer all three paths for the processing of scRNAseq data so it remains up to the user to decide which pipeline workflow is chosen for a particular analysis question.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
