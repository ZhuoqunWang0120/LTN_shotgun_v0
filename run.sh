#!/usr/bin/bash

export WORK_DIR=...

mkdir $WORK_DIR/cache
mkdir $WORK_DIR/results

# process data

$WORK_DIR/src/data/process_data_full.R $WORK_DIR 100
$WORK_DIR/src/data/process_data_simulation.R $WORK_DIR 100

# simulation: covariance estimation

$WORK_DIR/src/simulation/covariance_estimation/dtm_data_clrcov.R $WORK_DIR 1000000
sbatch $WORK_DIR/src/slurm/dtm_sim.q
sbatch $WORK_DIR/src/slurm/ln_sim1.q
sbatch $WORK_DIR/src/slurm/ln_sim2.q
sbatch $WORK_DIR/src/slurm/ln_sim3.q
sbatch $WORK_DIR/src/slurm/dtm_fit_parse.q
sbatch $WORK_DIR/src/slurm/ln_fit_parse.q
for i in `seq 1 100`;do
$WORK_DIR/src/simulation/covariance_estimation/fit_COAT.R $WORK_DIR $i
done
declare -a lams=(0.1 1 10)
for lambda in "${lams[@]}"; do
$WORK_DIR/src/simulation/covariance_estimation/collect_results_use.R $WORK_DIR $lambda 100 100
done
$WORK_DIR/src/simulation/covariance_estimation/clrcov_sensitivity_plot.R $WORK_DIR 100

# simulation: cross-group comparison 
sbatch $WORK_DIR/src/slurm/single_ltn_reffcov2.q
sbatch $WORK_DIR/src/slurm/multi_ltn_reffcov2.q
sbatch $WORK_DIR/src/slurm/single_dirfactor_r2.q 
sbatch $WORK_DIR/src/slurm/multi_dirfactor_r2.q 
$WORK_DIR/src/simulation/cross_group_comparison/roc.R $WORK_DIR 10
$WORK_DIR/src/simulation/cross_group_comparison/plot.R $WORK_DIR

# simulation: tree mis-specification
declare -a lams=(0.1 1 3 10 30 100)
${WORK_DIR}/src/slurm/mkall_tree_q.sh
${WORK_DIR}/src/slurm/mkall_map2T1.sh
$WORK_DIR/src/simulation/tree_misspecification/tree_robustness_replicated.R $WORK_DIR
$WORK_DIR/src/simulation/cross_group_comparison/plot.R $WORK_DIR

# case study
sbatch $WORK_DIR/src/slurm/application.q
$WORK_DIR/src/application/figures_tables.R $WORK_DIR


