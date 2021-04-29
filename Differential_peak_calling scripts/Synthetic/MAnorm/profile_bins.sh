#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "logs.err"
#SBATCH --output "logs.out"
#SBATCH --cpus-per-task=10
#SBATCH --partition="normal"
#SBATCH --mem=100000
#SBATCH --ntasks 1
#SBATCH --job-name "profile_bins_H3K4"

source /wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/pythonenv/python2.xenv/bin/activate

profile_bins --peaks=../../MACS2/Control_rep1/Control_rep1_DBpeaks.broadPeak,../../MACS2/Control_rep2/Control_rep2_DBpeaks.broadPeak,../../MACS2/Control_rep3/Control_rep3_DBpeaks.broadPeak,../../MACS2/Control_rep4/Control_rep4_DBpeaks.broadPeak,../../MACS2/Control_rep5/Control_rep5_DBpeaks.broadPeak,../../MACS2/Chip_rep1/Chip_rep1_DBpeaks.broadPeak,../../MACS2/Chip_rep2/Chip_rep2_DBpeaks.broadPeak,../../MACS2/Chip_rep3/Chip_rep3_DBpeaks.broadPeak,../../MACS2/Chip_rep4/Chip_rep4_DBpeaks.broadPeak,../../MACS2/Chip_rep5/Chip_rep5_DBpeaks.broadPeak \
--reads=../../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.1.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.2.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.3.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.4.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_control_merged_sorted.5.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.1.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.2.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.3.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.4.bam.bed,../../data/MCF7_H3K36me3_pooled_sorted_downsampled_merged_sorted.5.bam.bed \
--labs=s11,s21,s31,s41,s51,s12,s22,s32,s42,s52 -n Simulated
EOF
