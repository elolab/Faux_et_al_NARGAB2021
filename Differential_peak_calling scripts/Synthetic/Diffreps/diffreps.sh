#!/bin/sh

###############################################
# launch diffreps on SLURM                    #
###############################################
#                                             #
# #                                           #
#                                             #
###############################################



sbatch <<EOF
#!/bin/bash
#SBATCH --error "hermmanndiffreps.err"
#SBATCH --output "hermmanndiffreps.out"
#SBATCH --partition="normal"
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 6
#SBATCH --job-name diffReps

TREATMENT="../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.1.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.2.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.3.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.4.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.5.bam.bed" 
CONTROL="../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.1.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.2.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.3.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.4.bam.bed ../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.5.bam.bed"  
JOB_NAME="Diffreps"


/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/diffreps-master/bin/diffReps.pl -tr \${TREATMENT}  -co \${CONTROL}  -gn hg19 -re results/test_diffReps_simulated_nb.txt --pval=1

EOF
