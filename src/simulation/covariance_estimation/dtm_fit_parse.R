#!/usr/bin/env Rscript

library(GetoptLong,quietly = T)
library(mvtnorm,quietly = T)
niter=10000
nmc=10^6
ntop = 100
GetoptLong(
  "niter=i","number of Gibbs iterations.",
  "nmc=i","number of Monte Carlo samples.",
  "WORK_DIR=s","working directory.",
  "ntop=i","number of OTUs"
)
source(paste0(WORK_DIR,'/src/LTN/gibbs_hp_cov_mu.R'))
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
run = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
modelCov = 1
if (modelCov == 1){ # skip repetitive runs
  param = read.csv(paste0(WORK_DIR,'/config/params_cov_est_reorder.csv'))[run,]
  lambda = param$lambda
  SEED = param$seed
  i = SEED
  mu_var = param$mu_var
  mu_cov = mu_var * diag(ntop - 1)
  if (lambda==0){
    lambda='hp'
  }
  filenam = paste0('i',SEED,'_lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds')
  result_dir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/ltn/")
  system(paste0('mkdir -p ',result_dir))
  datadir=paste0(WORK_DIR,'/cache/covariance_estimation/dtm/')
  if (!file.exists(paste0(result_dir,filenam))){
    input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
    tree=input_data$tree
    dat_i=readRDS(paste0(datadir,'sim_ntop',ntop,'seed',i,'.rds'))
    yyl=dat_i[1:2]
    Y=yyl$Y
    YL=yyl$YL
    N=nrow(Y)
    p=ncol(Y)
    K=p+1
    S=1
    st1<-system.time(t1<-try(gibbs1<-gibbs_ltn_hp_cov_mu(niter=niter,YL=YL,Y=Y,lambda=lambda,hp_cov_mu=mu_cov,verbose=T)))
try(print(gibbs1$iteration))
# print(t1)
 if ("try-error" %in% class(t1)) {
	    saveRDS(gibbs1,paste0(result_dir,'/S/','gibbs_i',SEED,'_lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds'))
 stop('error in gibbs sampling') 
 }
} 
print(st1)
saveRDS(gibbs1,paste0(result_dir,'/S/','gibbs_i',SEED,'_lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds'))

     samp_seq=ceiling(niter/2):niter
     MU=gibbs1$MU
     OMEGA=gibbs1$OMEGA
     mu=apply(MU[,samp_seq],1,mean)
     omega=Reduce('+',OMEGA[samp_seq])/length(OMEGA[samp_seq])
     sigma=solve(omega)
     rm(gibbs1)
     st2<-system.time(clrcov_ltn<-clrcov_sim_log(mu,sigma,tree,nmc,F,NULL))
     if (sum(is.na(clrcov_ltn))+sum(is.infinite(clrcov_ltn))==0){
       saveRDS(clrcov_ltn,paste0(result_dir,filenam))
    } else {
      warning('NA or Inf in clrcov')
    }
  }
  

