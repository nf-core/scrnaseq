workflow STARSOLO {

        // build star index
    if (params.aligner == 'star') {

      if (!params.star_index && params.fasta){
        STAR_GENOMEGENERATE( genome_fasta_makeSTARindex, gtf_makeSTARindex )
      }


        
    }
}