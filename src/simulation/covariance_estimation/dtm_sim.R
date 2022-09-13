#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
ntop = as.numeric(argv[2])
i = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/")
filename = paste0('sim_ntop',ntop,'seed',i,'.rds')
if (!file.exists(paste0(cachedir,filename))){
  nsim=200
  total=10^5
  input_data1=readRDS(paste0(cachedir,"dtm_diab_MoM",ntop,".rds"))
  input_data2=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
  theta=as.vector(input_data1[[1]])
  tau=as.vector(input_data1[[2]])
  tree=input_data2$tree
  set.seed(i)
  sim_i=dtm_sim(nsim,tree,theta,tau,total)
  saveRDS(sim_i,paste0(cachedir,filename))
}