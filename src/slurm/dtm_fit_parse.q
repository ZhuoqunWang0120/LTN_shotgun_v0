#!/bin/bash
#SBATCH --array=1-1200
#SBATCH --nodes=1

#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dtm_fit_parse-%a.out
#SBATCH --error=logs/dtm_fit_parse-%a.err
srun $WORK_DIR/src/simulation/covariance_estimation/dtm_fit_parse.R --WORK_DIR $WORK_DIR --niter 10000 --ntop 100 