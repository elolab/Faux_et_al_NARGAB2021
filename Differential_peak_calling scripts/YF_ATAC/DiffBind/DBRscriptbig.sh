#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "DB.err"
#SBATCH --output "DB.out"
#SBATCH --cpus-per-task=5
#SBATCH --partition="long"
#SBATCH --mem=50000
#SBATCH --ntasks 1
#SBATCH --job-name "ATAC_DB"

module add R/3.6.1

Rscript --max-ppsize=50000 Diffbinds_backbone.R

EOF


