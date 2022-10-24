/*
 * Alignment with Cellranger Arc
 */

// include {CELLRANGER_ARC_MKGTF} from "../../modules/local/cellranger_arc/mkgtf/main.nf"
// include {CELLRANGER_ARC_MKREF} from "../../modules/nf-core/modules/cellranger/mkref/main.nf"
include {CELLRANGER_ARC_COUNT} from "../../modules/local/cellranger_arc/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ARC_ALIGN {
    take:
        fasta
        gtf
        cellranger_arc_index
        ch_folders

    main:
        ch_versions = Channel.empty()

        assert cellranger_arc_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--genome_fasta') and a gtf file ('--gtf')."

        if (!cellranger_arc_index) {
            /*
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_ARC_MKGTF( gtf )
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKGTF.out.versions)

            CELLRANGER_ARC_MKGTF.out.gtf.view()
            // Make reference genome
            CELLRANGER_ARC_MKREF( fasta, CELLRANGER_ARC_MKGTF.out.gtf, "cellranger_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_ARC_MKREF.out.versions)
            cellranger_arc_index = CELLRANGER_ARC_MKREF.out.reference
            */
        }

        // Obtain read counts
        /*CELLRANGER_ARC_COUNT (
            ch_folders,
            cellranger_arc_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_ARC_COUNT.out.versions)
        */
    emit:
        ch_versions
        cellranger_arc_out  = CELLRANGER_ARC_COUNT.out.outs
}


def create_lib_csv(){
    echo "fastqs,sample,library_type" > $ANA_PATH/$sample/lib.csv
    fastq_path=$DATA_PATH/${sample}_ATAC
    fastq_heads=`ls $fastq_path | grep -oP '.+(?=_S[0-9]+_L001_R1_001.fastq.gz)'`
    for fastq_head in $fastq_heads; do
        echo "$fastq_path,$fastq_head,Chromatin Accessibility" >> $ANA_PATH/$sample/lib.csv
    done
    fastq_path=$DATA_PATH/${sample}_GEX
    fastq_head=`ls $fastq_path | grep -oP '.+(?=_S[0-9]+_L001_R1_001.fastq.gz)'`
    echo "$fastq_path,$fastq_head,Gene Expression" >> $ANA_PATH/$sample/lib.csv
}