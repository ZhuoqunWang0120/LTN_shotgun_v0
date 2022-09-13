#!/usr/bin/bash

lam=$1
echo \#\!/bin/bash
echo \#SBATCH --array=1-100
echo \#SBATCH --nodes=1
echo \#SBATCH --mem=4G
echo \#SBATCH --ntasks=1
echo \#SBATCH --cpus-per-task=1
echo \#SBATCH --output=logs/balanced$lam-%a.out
echo \#SBATCH --error=logs/balanced$lam-%a.err
echo srun $WORK_DIR/R/tree/tree_robustness_ex0_balanced.R $WORK_DIR $lam
