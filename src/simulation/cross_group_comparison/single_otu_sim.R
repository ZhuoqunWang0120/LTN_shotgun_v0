#!/usr/bin/env Rscript

# library(LTN)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
i = ceiling(id / 2)
h = 2 * i - id
# i=as.numeric(argv[2])
# h=as.numeric(argv[3])
ntop = 100
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
datadir=paste0(WORK_DIR,"/cache/cross_group_comparison/single_otu/")
system(paste0('mkdir -p ',datadir))
filenam=paste0('sim',i,'H',h,'.rds')
if (!file.exists(paste0(datadir,filenam))){
  input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
  cnt=input_data$cnt
  tree=input_data$tree
  K=ncol(cnt)
  N=nrow(cnt)
  set.seed(i)
  group2=sample(N,ceiling(N/2),replace = F)
  otu20=input_data$top20_ra
  if (h==1){
    otu_modify=sample(otu20,1)
    cnt[group2,otu_modify]=3*cnt[group2,otu_modify]
  }
  Xtest=matrix(0,N,1)
  Xtest[group2,]=1
  yyl=seqtab2y(cnt,tree)
  saveRDS(list(Y=yyl$Y,YL=yyl$YL,cnt=cnt,Xtest=Xtest,group2=group2),paste0(datadir,filenam))
}

