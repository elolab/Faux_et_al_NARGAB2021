#!/bin/sh

###############################################
# launch PePr on SLURM                        #
###############################################



sbatch <<EOF
#!/bin/bash
#SBATCH --error "peprno.err"
#SBATCH --output "peprno.out"
#SBATCH --partition="normal"
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --job-name pepr
source /wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/pythonenv/python2.xenv/bin/activate

TREATMENT="../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.1.bam,../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.2.bam,../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.3.bam,../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.4.bam,../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.5.bam" 
TREATMENT2="../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.1.bam,../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.2.bam,../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.3.bam,../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.4.bam,../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.5.bam"  
JOB_NAME="test_Simulated"


PePr -c \${TREATMENT} --chip2 \${TREATMENT2} -n \${JOB_NAME} -f bam --normalization no --peaktype broad --diff --output-directory Without_normalization --threshold=0.5

EOF
