# nf-core/scrnaseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## Kallisto & Bustools Results

See [Kallisto](https://pachterlab.github.io/kallisto/about) for details about Kallisto and [Bustools](https://bustools.github.io/) for more information on BusTools. 

The pipeline can analyze data from single cell rnaseq experiments and generates a set of folders with respective outputs from various steps of the analysis. For a detailed summary what the pipeline does specifically, please follow the [excellent tutorial](https://www.kallistobus.tools/getting_started.html) that also describes specific steps for downstream analysis of the generated matrices. 

**Output directory: `results/kallisto`**

* `kallisto_gene_map`
  * Contains the converted GTF gene map that is used by BusTools for downstream analysis
* `kallisto_index`
  * Contains the index of the supplied (genome/transcriptome) FastA file
* `bustools_counts`
  * Contains two subdirectories
    * `eqcount`: Containing the Transcript Compatibility Count (TCC) Matrix (`tcc.mtx`)
    * `genecount`: Containing the Gene Count Matrix (`gene.mtx`)

For details on how to load these into R and perform further downstream analysis, please refer to the [BusTools HowTo](https://github.com/BUStools/getting_started/blob/master/getting_started.ipynb). 



## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
