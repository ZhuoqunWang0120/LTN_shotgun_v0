#!/bin/bash
#SBATCH --array=1-100
#SBATCH --nodes=1

#SBATCH --mem=1G       
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/ln_sim1-%a.out
#SBATCH --error=logs/ln_sim1-%a.err
srun $WORK_DIR/src/simulation/covariance_estimation/ln_sim.R $WORK_DIR 1 100 