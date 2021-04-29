#!/bin/sh

sbatch <<EOF
#!/bin/bash
#SBATCH --error "Rlogs.err"
#SBATCH --output "Rlogs.out"
#SBATCH --cpus-per-task=4
#SBATCH --partition="normal"
#SBATCH --ntasks 1
#SBATCH --job-name "MA_normATAC"

module add R

Rscript --max-ppsize=100000 MAnorm2.r
EOF
