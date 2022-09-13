#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
nsim1=argv[3]
ntop = argv[4]
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
if (lambda=='0'){lambda='hp'}
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/")
system(paste0('mkdir -p ',result_dir))

cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/ltn/")
for (mu_var in c(1,5,50,100)){
  risk=matrix(-1,4,4,dimnames = list(c("Frobenius","L1","L_inf","spectral"),c("LN-hub","LN-block","LN-sparse","DTM")))
  for (m in c('hub','block','sparse')){
    cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/ltn/")
    files=grep(paste0(m,'_lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds'),list.files(cachedir,full.names = T),value = T)
    loss=do.call(rbind,lapply(1:nsim1, function(x){
      clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/ln/",m,"clrcov_",x,'ntop',ntop,".rds"))
      clrcov1=readRDS(grep(paste0('i',x,'_'),files,value=T))
      matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
    risk[,paste0("LN-",m)]=(colSums(loss)/nrow(loss))
  }
  cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/ltn/")
  files=grep(paste0('lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds'),list.files(cachedir,full.names = T),value = T)
  clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/dtm/clrcov_true",ntop,".rds"))
  loss=do.call(rbind,lapply(1:nsim1, function(x){
    clrcov1=readRDS(grep(paste0('i',x,'_'),files,value=T))
    matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
  risk[,"DTM"]=(colSums(loss)/nrow(loss))
  saveRDS(risk,paste0(result_dir,"/loss_LTN_lambda_",lambda,'_mu_var',mu_var,".rds"))
}


loss = list()
for (mu_var in c(1,5,50,100)){
  loss[[paste0('mu_var',mu_var)]] = readRDS(paste0(result_dir,"/loss_LTN_lambda_",lambda,'_mu_var',mu_var,".rds"))
}

print(lapply(loss, function(x)round(x,2)))
loss = list()
for (mu_var in c(1,5,50,100)){
  loss[[paste0('mu_var',mu_var)]] = readRDS(paste0(result_dir,"/loss_LTN_lambda_",lambda,'_mu_var',mu_var,".rds"))
}
# COAT
risk=matrix(-1,4,4,dimnames = list(c("Frobenius","L1","L_inf","spectral"),c("LN-hub","LN-block","LN-sparse","DTM")))
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/coat/")
clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/dtm/clrcov_true100.rds"))
files=list.files(cachedir,full.names = T)
loss=do.call(rbind,lapply(files, function(x){matloss(cov2cor(readRDS(x))-cov2cor(clrcov0))}))
risk[,'DTM']=(colSums(loss)/nrow(loss))
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/coat/")
for (m in c('hub','block','sparse')){
  files=grep(m,list.files(cachedir,full.names = T),value = T)
  loss=do.call(rbind,lapply(1:nsim1, function(x){
    clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/ln/",m,"clrcov_",x,'ntop',ntop,".rds"))
    clrcov1=readRDS(grep(paste0('est',x,'ntop',ntop,'.rds'),files,value=T))
    matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
  risk[,paste0("LN-",m)]=(colSums(loss)/nrow(loss))
}
saveRDS(risk,paste0(result_dir,"/loss_COAT_ntop",ntop,".rds"))

print(round(risk,2))