#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --nodes=1

#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/m_ltn_2-%a.out
#SBATCH --error=logs/m_ltn_2-%a.err

reffcov=2
lambda=10
srun $WORK_DIR/src/simulation/cross_group_comparison/multi_otu_sim.R $WORK_DIR logs/m_sim$SLURM_ARRAY_TASK_ID
srun $WORK_DIR/src/simulation/cross_group_comparison/fit.R --reffcov $reffcov --lambda $lambda --scenario multi_otu --niter 10000 --WORK_DIR $WORK_DIR logs/multi_ltn_$SLURM_ARRAY_TASK_ID

