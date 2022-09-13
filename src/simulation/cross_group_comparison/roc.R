#!/usr/bin/env Rscript

library(ROCR)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
for (s in c('single_otu','multi_otu')){
  system(paste0('mkdir -p ',WORK_DIR,'/results/cross_group_comparison/',s))
  for (m in c('sparse','diagonal')){
    cachedir=paste0(WORK_DIR,"/cache/cross_group_comparison/",s,"/LTN/pjap/")
    files0=grep(paste0('H0_',m,'_lambda',lambda,'.rds'),list.files(cachedir,full.names = T),value = T)
    files1=grep(paste0('H1_',m,'_lambda',lambda,'.rds'),list.files(cachedir,full.names = T),value = T)
    pjaps=c(sapply(files0, readRDS),sapply(files1, readRDS))
    labels=c(rep(0,length(files0)),rep(1,length(files1)))
    if (length(unique(labels))>1){
      if (lambda=='0' & m=='sparse'){lambda1='GammaPrior'}else{if (m=='diagonal'){lambda1='NA'}else{lambda1=lambda}}
      roc=performance(prediction(pjaps,labels),'tpr','fpr')
      saveRDS(roc,paste0(WORK_DIR,'/results/cross_group_comparison/',s,'/LTN_',m,'_lambda',lambda1,'.rds'))
    }
  }
}
