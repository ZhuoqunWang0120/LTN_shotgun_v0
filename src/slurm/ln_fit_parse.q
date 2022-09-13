#!/bin/bash
#SBATCH --array=1201-3600
#SBATCH --nodes=1

#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/ln_fit_parse-%a.out
#SBATCH --error=logs/ln_fit_parse-%a.err
srun $WORK_DIR/src/simulation/covariance_estimation/ln_fit_parse.R --WORK_DIR $WORK_DIR --nmc 1000000 --niter 10000 --ntop 100 logs/ln_fit_parse_$SLURM_ARRAY_TASK_ID

