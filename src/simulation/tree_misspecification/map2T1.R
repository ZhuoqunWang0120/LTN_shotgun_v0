#!/usr/bin/env Rscript
argv = commandArgs(TRUE)
WORK_DIR = argv[1]
lambda = as.numeric(argv[2])
tree_structure = as.numeric(argv[3])
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(mvtnorm)
set.seed(1)
SRC_DIR = paste0(WORK_DIR,'/src/simulation/tree_misspecification/')
RESULT_DIR = paste0(WORK_DIR,'/cache/tree/')
source(paste0(SRC_DIR,'tree_structure.R'))
system(paste0('mkdir ',RESULT_DIR,'map2T1/'))
K = 64
nmc = 10^6
if (tree_structure == 1){
  # right2left
if (!(file.exists(paste0(RESULT_DIR,'map2T1/est_right2left_','lambda',lambda,'K',K,'_',id,'.rds')))){
  estpara = readRDS(paste0(RESULT_DIR,'est_para_right_lambda',lambda,'K',K,'_',id,'.rds'))
  tictoc::tic()
  para_original = right2left(estpara$mu_right,estpara$sigma_right,nmc)
  tictoc::toc()
  saveRDS(para_original, paste0(RESULT_DIR,'map2T1/est_right2left_','lambda',lambda,'K',K,'_',id,'.rds'))
}}
if (tree_structure == 2){
  # balanced2left
  if (!(file.exists(paste0(RESULT_DIR,'map2T1/est_balanced2left_','lambda',lambda,'K',K,'_',id,'.rds')))){
	estpara = readRDS(paste0(RESULT_DIR,'est_para_balanced_lambda',lambda,'K',K,'_',id,'.rds'))
  tictoc::tic()
  para_original = balanced2left(estpara$mu_balanced,estpara$sigma_balanced,nmc)
  tictoc::toc()
  saveRDS(para_original, paste0(RESULT_DIR,'map2T1/est_balanced2left_','lambda',lambda,'K',K,'_',id,'.rds'))
}
}






