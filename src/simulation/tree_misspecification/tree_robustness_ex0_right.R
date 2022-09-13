#!/usr/bin/env Rscript
argv = commandArgs(TRUE)
WORK_DIR = argv[1]
lambda = as.numeric(argv[2])
system(paste0('mkdir ',WORK_DIR,'/cache/tree/'))
cachedir = paste0(WORK_DIR,'/cache/tree/')
source(paste0(WORK_DIR,'/src/simulation/tree_misspecification/tree_structure.R'))
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
source(paste0(WORK_DIR,"/src/LTN/gibbs_glasso.R"))
source(paste0(WORK_DIR,"/src/LTN/gibbs_hp_cov_mu.R"))

library(mvtnorm)
# true tree is left-skewed
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(id)
K = 64
n = 200
mu = rep(2,K-1)
sigma = diag(K-1)
psi = rmvnorm(n = n,mean = mu,sigma = sigma)
nodevals = theta2nodevals_special(plogis(psi),'left')
prob = node2leaf_special(nodevals, 'left')
if (sum(abs(rowSums(prob)-1)>10^(-10))>0){
  stop('incorrect leaf prob: do not sum to 1')
}
cnt = t(apply(prob, 1, function(x) rmultinom(1,10^5,x)))
saveRDS(cnt, paste0(cachedir,'cnt_left_lambda',lambda,'K',K,'_',id,'.rds'))
# fit LTN with right-skewed tree
niter = 10^4
yyl = leaf2node_special(cnt, 'right')
if (!file.exists(paste0(cachedir,'gibbs_right_lambda',lambda,'K',K,'_',id,'.rds'))){

gibbs1 = gibbs_ltn_hp_cov_mu(niter=niter,YL=yyl$leftval,Y=yyl$nodeval,SEED=1,lambda=lambda)
saveRDS(gibbs1, paste0(cachedir,'gibbs_right_lambda',lambda,'K',K,'_',id,'.rds'))
samp_seq=(niter/2+1):niter
MU=gibbs1$MU
OMEGA=gibbs1$OMEGA
rm(gibbs1)
mu=apply(MU[,samp_seq],1,mean)
omega=Reduce('+',OMEGA[samp_seq])/length(OMEGA[samp_seq])
sigma=solve(omega)
saveRDS(list(mu_right=mu,sigma_right=sigma), paste0(cachedir,'est_para_right_lambda',lambda,'K',K,'_',id,'.rds'))
}


