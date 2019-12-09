# ![nf-core/scrnaseq](docs/images/nfcore-scrnaseq_logo.png)

**A fully automated Nextflow pipeline for Droplet-based (e.g. 10x Genomics) single-cell RNA-Seq data**.

[![Build Status](https://github.com/nf-core/scrnaseq/workflows/scrnaseq%20CI/badge.svg)](https://github.com/nf-core/scrnaseq/workflows/scrnaseq%20CI/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/scrnaseq.svg)](https://hub.docker.com/r/nfcore/scrnaseq)

[![Join us on Slack](https://img.shields.io/badge/slack-nfcore/scrnaseq-blue.svg)](https://nfcore.slack.com/channels/scrnaseq)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

Work in progress - this is a community effort in building a pipeline capable to support:

* Alevin + AlevinQC
* STARSolo
* Kallisto + BUStools

## Documentation

The nf-core/scrnaseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

The `nf-core/scrnaseq` was initiated by [Peter J. Bailey](https://github.com/PeterBailey) (Salmon Alevin, AlevinQC) with major contributions from [Olga Botvinnik](https://github.com/olgabot) (STARsolo, Testdata) and [Alex Peltzer](https://github.com/apeltzer) (Kallisto/BusTools workflow).

## Citation

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).

The basic benchmarks that were used as motivation for incorporating the three available modular workflows can be found in [this publication](https://www.biorxiv.org/content/10.1101/673285v2).

We offer all three paths for the processing of scRNAseq data so it remains up to the user to decide which pipeline workflow is chosen for a particular analysis question.
