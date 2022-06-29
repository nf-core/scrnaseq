//
// Check input samplesheet and get read channels
//

//TODO --> add skip_fastqc to params

include { FASTQC } from '../../modules/nf-core/modules/fastqc/main'

workflow FASTQC_CHECK {
  take:
  ch_fastq

  main:
  ch_fastq
      .map { ch -> [ ch[0], ch[1] ] }
      .set { ch_fastq }

  /*
   * FastQ QC using FASTQC
   */
  fastqc_zip     = Channel.empty()
  fastqc_html    = Channel.empty()
  fastqc_multiqc = Channel.empty()
  fastqc_version = Channel.empty()
  
  FASTQC ( ch_fastq )
  fastqc_zip     = FASTQC.out.zip
  fastqc_html    = FASTQC.out.html

  fastqc_zip
      .map { it -> [ it[1] ] }
      .set { fastqc_zip_only }
  fastqc_html
      .map { it -> [ it[1] ] }
      .set { fastqc_html_only }

  fastqc_multiqc = fastqc_multiqc.mix( fastqc_zip_only, fastqc_html_only )
  fastqc_version = FASTQC.out.versions

  emit:
  fastqc_zip
  fastqc_html
  fastqc_version
  fastqc_multiqc
}
