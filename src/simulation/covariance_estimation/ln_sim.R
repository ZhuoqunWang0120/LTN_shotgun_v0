#!/usr/bin/env Rscript

library(philr,quietly = T)
library(mvtnorm,quietly = T)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
modelCov=as.numeric(argv[2])
ntop = as.numeric(argv[3])
SEED = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source(paste0(WORK_DIR,"/src/simulation/covariance_estimation/ln_graph.R"))
muvar=16
modelcovs=c('hub','block','sparse')
modelcov=modelcovs[modelCov]
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
input_data=readRDS(paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))
tree=input_data$tree
datadir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/")
system(paste0('mkdir -p ',datadir))
set.seed(SEED)
parasimu=generateDataCell(n=200,modelCov=modelCov,p1=(ntop - 1),p2=100,p3=200,nRep=1)
mu=rnorm((ntop - 1),0,sqrt(muvar))
saveRDS(list(mu=mu,sigma=solve(parasimu$sigma1)),paste0(datadir,modelcov,'para_',SEED,'ntop',ntop,'.rds'))
# use sigma1 as omega
eta=rmvnorm(200,mu,solve(parasimu$sigma1))
colnames(eta)=paste0('n',1:(ntop - 1))
tree$node.label=paste0('n',1:(ntop - 1))
philinv=philrInv(df.ilrp = eta,tree=tree,V=NULL)
cnt=t(apply(philinv,1,function(x){rmultinom(1,10^5,x)}))
colnames(cnt)=tree$tip.label
saveRDS(cnt,paste0(datadir,modelcov,'cnt_',SEED,'ntop',ntop,'.rds'))
yyl=seqtab2y(cnt,tree)
Y=yyl$Y
YL=yyl$YL
saveRDS(list(Y=Y,YL=YL),paste0(datadir,modelcov,'yyl_',SEED,'ntop',ntop,'.rds'))

# clrcov via transformation
sig=solve(parasimu$sigma1)
K=ntop
p=K-1
ilrE=diag(p)
colnames(ilrE)=tree$node.label
E=philrInv(ilrE,tree)
H=t(apply(E, 1, function(x){clrp(x,rep(1,K))})) # clrp=clr when p=1
G=t(H)%*%sig%*%H
saveRDS(G,paste0(datadir,modelcov,'clrcov_',SEED,'ntop',ntop,'.rds'))