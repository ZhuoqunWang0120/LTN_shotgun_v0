#!/bin/bash
#SBATCH --array=1-100
#SBATCH --nodes=1

#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dtm_sim-%a.out
#SBATCH --error=logs/dtm_sim-%a.err
srun $WORK_DIR/src/simulation/covariance_estimation/dtm_sim.R $WORK_DIR 100 
