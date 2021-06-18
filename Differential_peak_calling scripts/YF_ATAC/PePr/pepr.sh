#!/bin/sh

###############################################
# launch diffreps on SLURM                    #
###############################################



sbatch <<EOF
#!/bin/bash
#SBATCH --error "peprH3K36me3.err"
#SBATCH --output "peprH3K36me3.out"
#SBATCH --partition="long"
#SBATCH --time=4-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --job-name peprATAC

source /wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/pythonenv/python2.xenv/bin/activate

TREATMENT="../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839476_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839477_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839478_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839479_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839480_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839481_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839482_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839483_marked_woblr.bam" 
CONDITION="../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839491_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839492_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839493_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839494_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839495_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839496_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839497_marked_woblr.bam,../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839498_marked_woblr.bam"  
JOB_NAME="ATAC_P"



PePr -c \${TREATMENT} --chip2 \${CONDITION} -n \${JOB_NAME} -f bam --peaktype=broad --diff --output-directory results  --threshold 0.1

EOF
