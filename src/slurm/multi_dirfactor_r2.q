#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --nodes=1

#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/m_dirfactor_r2-%a.out
#SBATCH --error=logs/m_dirfactor_r2-%a.err

r=2
srun $WORK_DIR/src/simulation/cross_group_comparison/dirfactor/R/fit_multi_otu.R $WORK_DIR $r 100000 logs/m_dirfactor_r${r}_$SLURM_ARRAY_TASK_ID

