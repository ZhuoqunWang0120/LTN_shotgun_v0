#!/usr/bin/env Rscript
library(GetoptLong)
ntop=100
niter=10000
adjust=T
gm=100
pnull=0.5
r=0
save_alpha_only=T
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
i = ceiling(id / 2)
h = 2 * i - id
GetoptLong(
  "reffcov=i","random effect covariance, 1: diagonal, 2: sparse.",
  "h=i","null or alternative.",
  "lambda=f","lambda, 0 represents Gamma prior on it.",
  "i=i","seed of simulation.",
  "scenario=s","simulation scenario.",
  "niter=i","number of Gibbs iterations.",
  "WORK_DIR=s","working directory."
)


# library(LTN)
source(paste0(WORK_DIR, '/src/LTN/utils.R'))
source(paste0(WORK_DIR, '/src/LTN/mixed_effects.R'))
cachedir=paste0(WORK_DIR,"/cache/cross_group_comparison/",scenario,"/LTN/")
model_index=reffcov
gprior_m=gm
gibbsdir=paste0(cachedir,'/gibbs/')
pmapdir=paste0(cachedir,'/pmap/')
pjapdir=paste0(cachedir,'/pjap/')
timedir=paste0(cachedir,'/time/')
sdir=paste0(cachedir,'/SEED/')
randmodel=c('norand','diagonal','sparse')[model_index+1]
try(system(paste0('mkdir -p ',cachedir)))
try(system(paste0('mkdir -p ',pmapdir)))
try(system(paste0('mkdir -p ',pjapdir)))
try(system(paste0('mkdir -p ',timedir)))
try(system(paste0('mkdir -p ',sdir)))
filenam=paste0('i',i,'H',h,'_',randmodel,'_lambda',lambda,'.rds')
if (!file.exists(paste0(cachedir,'/pjap/',filenam))){
  input_data1=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
  input_data2=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/",scenario,"/sim",i,"H",h,".rds"))
  Y=input_data2$Y
  YL=input_data2$YL
  N=nrow(Y)
  p=ncol(Y)
  K=p+1
  g=input_data1$g
  Xtest=input_data2$Xtest
  Xadjust=input_data1$Xadjust
  grouplabel=input_data1$grouplabel
  SEED=1
  st<-system.time(t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,niter=niter,adjust=adjust,reff=T,reffcov=reffcov,SEED=SEED,save_alpha_only=save_alpha_only,gprior_m=gm,pnull=pnull,a1a2='none',lambda_fixed=lambda,verbose=T)))
  if ("try-error" %in% class(t)) {
    # SEED=SEED+1
    warning('numerical issues')
    saveRDS(list(SEED,t,id),paste0(sdir,'/',filenam))
    # t<-try(gibbs<-gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,niter=niter,adjust=adjust,reff=T,reffcov=reffcov,SEED=SEED,save_alpha_only=save_alpha_only,gprior_m=gm,pnull=pnull,a1a2='none',lambda_fixed=lambda))
  }
  BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
  BETAMAT1=do.call(rbind,BETA1)
  PMAP=matrix(apply(BETAMAT1[(niter/2+1):niter,],2,function(x){sum(x!=0)})/(niter/2),nrow = 1)
  saveRDS(PMAP,paste0(pmapdir,'/',filenam))
  all0=rowSums(BETAMAT1[(niter/2+1):niter,]!=0)
  PJAP=1-sum(all0==0)/(niter/2)
  saveRDS(PJAP,paste0(pjapdir,'/',filenam))
  saveRDS(st,paste0(timedir,'/',filenam))
  print(st)
}
