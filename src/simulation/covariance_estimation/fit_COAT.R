#!/usr/bin/env Rscript
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
SEED=as.numeric(argv[2])
i=SEED
source(paste0(WORK_DIR,"/src/simulation/covariance_estimation/COAT.R"))
#----------------------dtm------------------------------------
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/")
resultdir=paste0(cachedir,'coat/')
system(paste0('mkdir -p ',resultdir))
cnt=readRDS(paste0(cachedir,'sim',SEED,'.rds'))$cnt
cnt[cnt==0]=0.5
coat_est=coat(cnt)$sigma
saveRDS(coat_est,paste0(resultdir,'est',SEED,'.rds'))
#----------------------ilr------------------------------------
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/")
resultdir=paste0(cachedir,'coat/')
system(paste0('mkdir -p ',resultdir))
modelcovs=c('hub','block','sparse')
for (modelCov in 1:3){
  modelcov=modelcovs[modelCov]
  clrtrue=readRDS(paste0(cachedir,modelcov,'clrcov_',i,'.rds'))
  cnt=readRDS(paste0(cachedir,modelcov,'cnt_',i,'.rds'))
  cnt[cnt==0]=0.5
  coat_est=coat(cnt)$sigma
  saveRDS(coat_est,paste0(resultdir,modelcov,'_est',i,'.rds'))
  }