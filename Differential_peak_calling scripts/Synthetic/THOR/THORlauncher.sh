#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "THOR.err"
#SBATCH --output "THOR.out"
#SBATCH --partition="long"
#SBATCH --mem=20000
#SBATCH --ntasks 1
#SBATCH --job-name "THOR"

source /wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/pythonenv/python2.xenv/bin/activate

rgt-THOR --scaling-factors=1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0 --pvalue=1 --name="test_THOR-synthetic" THOR.config
EOF
