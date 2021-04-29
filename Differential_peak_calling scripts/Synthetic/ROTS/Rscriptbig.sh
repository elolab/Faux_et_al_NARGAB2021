#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "logs/base.err"
#SBATCH --output "logs/base.out"
#SBATCH --cpus-per-task=10
#SBATCH --partition="long"
#SBATCH --mem=100000
#SBATCH --ntasks 1
#SBATCH --job-name "RO_synth"

module add R/

Rscript --max-ppsize=100000 ROTS_backbone_base_1.0.5.R
EOF
