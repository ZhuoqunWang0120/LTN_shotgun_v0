#!/usr/bin/env Rscript

library(ROCR)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
for (scenario in c('single_otu','multi_otu')){
  datadir=paste0(WORK_DIR,'/cache/cross_group_comparison/',scenario,'/dirfactor/teststat/')
  resultdir=paste0(WORK_DIR,'/results/cross_group_comparison/',scenario,'/')
  try(system(paste0('mkdir -p ',resultdir)))
  for (r in c(2,5)){
    files=grep(paste0('r',r,'iter'),list.files(datadir),value=T)
    files=grep('norm',files,value=T)
    files0=grep('H0',files,value=T)
    files1=grep('H1',files,value=T)
    n0=length(files0)
    n1=length(files1)
    v0=sapply(files0,function(x){readRDS(paste0(datadir,x))})
    v1=sapply(files1,function(x){readRDS(paste0(datadir,x))})
    v=c(v0,v1)
    labels=c(rep(0,n0),rep(1,n1))
    #print(c(r,length(labels)))
    #saveRDS(list(v=v,labels=labels),paste0(datadir,'/summaries/v_labels_r',r,'.RData'))
    roc=performance(prediction(v,labels),'tpr','fpr')
    saveRDS(roc,paste0(resultdir,'/dirfactor_r',r,'.RData'))
  }
}

