#!/bin/bash
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/application-%a.out
#SBATCH --error=logs/application-%a.err
lambda=10
srun $WORK_DIR/src/application/application100.R $WORK_DIR 10000 $lambda


