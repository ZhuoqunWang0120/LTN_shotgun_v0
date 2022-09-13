#!/bin/bash
#SBATCH --array=1-100
#SBATCH --nodes=1

#SBATCH --mem=1G  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/ln_sim3-%a.out
#SBATCH --error=logs/ln_sim3-%a.err
srun $WORK_DIR/src/simulation/covariance_estimation/ln_sim.R $WORK_DIR 3 100 