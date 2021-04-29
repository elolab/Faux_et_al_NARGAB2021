#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "THOR.err"
#SBATCH --output "THOR.out"
#SBATCH --partition="long"
#SBATCH --time=48:00:00
#SBATCH --mem=20000
#SBATCH --ntasks 1
#SBATCH --job-name "THORATAC"

source /wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/pythonenv/python2.xenv/bin/activate

rgt-THOR --pvalue=0.05 --name="THOR-YF" THOR.config
EOF
