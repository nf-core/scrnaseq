workflow KALLISTO_BUSTOOLS {
        // kallisto
    if (params.aligner == 'kallisto') {

      // build index
      if (!params.kallisto_index){
        KALLISTO_INDEX( transcriptome_fasta_kallisto )
      }

      // TODO: kallisto genemap


    }
}