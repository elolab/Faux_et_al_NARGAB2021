#!/bin/sh

###############################################
# launch diffreps on SLURM                    #
###############################################



sbatch <<EOF
#!/bin/bash
#SBATCH --error "diffreps.err"
#SBATCH --output "diffreps.out"
#SBATCH --partition="long"
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --job-name diffReps



TREATMENT="../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839476_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839477_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839478_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839479_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839480_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839481_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839482_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839483_marked_woblr.bam.bed" 
CONTROL="../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839491_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839492_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839493_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839494_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839495_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839496_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839497_marked_woblr.bam.bed ../../../ChIP-seq_analysis/ATAC-seq/6-marked-sorted_BAM/SRR5839498_marked_woblr.bam.bed"  
JOB_NAME="ATAC_DR"

/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/diffreps-master/bin/diffReps.pl -tr \${TREATMENT} -co \${CONTROL} -gn hg19 -re results/diffReps_simulated_nb.txt --pval 0.1

EOF
