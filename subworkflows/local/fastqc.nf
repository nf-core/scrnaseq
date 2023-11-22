//
// Check input samplesheet and get read channels
//
include { FASTQC } from '../../modules/nf-core/fastqc/main'

workflow FASTQC_CHECK {
  take:
  ch_fastq

  main:

  if( params.aligner == 'cellrangerarc'){
    ch_fastq
        .map { ch -> [ ch[0], ch[3] ] }
        .set { ch_fastq }
  } else {
    ch_fastq
        .map { ch -> [ ch[0], ch[1] ] }
        .set { ch_fastq }
  }

  /*
   * FastQ QC using FASTQC
   */
  FASTQC ( ch_fastq )
  fastqc_zip     = FASTQC.out.zip
  fastqc_html    = FASTQC.out.html

  fastqc_zip
      .map { it -> [ it[1] ] }
      .set { fastqc_zip_only }
  fastqc_html
      .map { it -> [ it[1] ] }
      .set { fastqc_html_only }

  fastqc_multiqc = Channel.empty()
  fastqc_multiqc = fastqc_multiqc.mix( fastqc_zip_only, fastqc_html_only )
  fastqc_version = FASTQC.out.versions

  emit:
  fastqc_zip
  fastqc_html
  fastqc_version
  fastqc_multiqc
}
