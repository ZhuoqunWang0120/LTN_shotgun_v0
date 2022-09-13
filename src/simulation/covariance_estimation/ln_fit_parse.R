#!/usr/bin/env Rscript
# pkgs=c('GetoptLong','philr','mvtnorm','LTN')
pkgs = c('GetoptLong','philr','mvtnorm')
lapply(pkgs, require,character.only = TRUE)
niter=10000
nmc=10^6
ntop=100
GetoptLong(
  "niter=i","number of Gibbs iterations.",
  "nmc=i","number of Monte Carlo samples.",
  "WORK_DIR=s","working directory",
  "ntop=i","number of OTUs"
)
source(paste0(WORK_DIR,'/src/LTN/gibbs_hp_cov_mu.R'))
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
run = ceiling(id/3)
modelCov = id - (run - 1)*3
param = read.csv(paste0(WORK_DIR,'/config/params_cov_est1.csv'))[run,]
lambda = param$lambda
SEED = param$seed
mu_var = param$mu_var
mu_cov = mu_var * diag(ntop - 1)
if (lambda==0){
  lambda='hp'
}
modelcovs=c('hub','block','sparse')
modelcov=modelcovs[modelCov]
#source(paste0(WORK_DIR,"/src/experiments/covariance_estimation/gibbs.R"))
input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
tree=input_data$tree
datadir=paste0(WORK_DIR,'/cache/covariance_estimation/ln/')
result_dir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/ltn/")
system(paste0('mkdir -p ',result_dir))
filenam=paste0('i',SEED,'_',modelcov,'_lambda',lambda,'ntop',ntop,'mu_var',mu_var,'.rds')
yyl=readRDS(paste0(datadir,modelcov,'yyl_',SEED,'ntop',ntop,'.rds'))
Y=yyl$Y
YL=yyl$YL
N=nrow(Y)
p=ncol(Y)
K=p+1
S=1
if (!file.exists(paste0(result_dir,filenam))){
  st1<-system.time(t<-try(gibbs1<-gibbs_ltn_hp_cov_mu(niter=niter,YL=YL,Y=Y,lambda=lambda,hp_cov_mu=mu_cov,verbose=T)))
  while("try-error" %in% class(t)) {
    S=S+1
    warning('numerical issues')
    t<-try(gibbs1<-gibbs_ltn_hp_cov_mu(niter=niter,YL=YL,Y=Y,SEED=S,lambda=lambda,hp_cov_mu=mu_cov,verbose=T))
  }
  samp_seq=(niter/2+1):niter
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
