#!/usr/bin/bash


lam=$1
tree=$2
echo \#\!/bin/bash
echo \#SBATCH --array=1-100
echo \#SBATCH --nodes=1
echo \#SBATCH --mem=8G
echo \#SBATCH --ntasks=1
echo \#SBATCH --cpus-per-task=1
echo \#SBATCH --output=logs/map$lam-$tree-%a.out
echo \#SBATCH --error=logs/map$lam-$tree-%a.err
echo srun $WORK_DIR/R/tree/map2T1.R $WORK_DIR $lam $tree

