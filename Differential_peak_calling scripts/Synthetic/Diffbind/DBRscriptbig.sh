#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "DB.err"
#SBATCH --output "DB.out"
#SBATCH --cpus-per-task=10
#SBATCH --partition="normal"
#SBATCH --mem=50000
#SBATCH --ntasks 1
#SBATCH --job-name "Herrmann"

module add R

Rscript --max-ppsize=50000 Diffbinds_backbone.R
EOF
