#!/usr/bin/env Rscript
argv = commandArgs(TRUE)
WORK_DIR = argv[1]
lambda = as.numeric(argv[2])
cachedir = paste0(WORK_DIR,'/cache/tree/')
source(paste0(WORK_DIR,'/src/simulation/tree_misspecification/tree_structure.R'))
source(paste0(WORK_DIR,"/src/LTN/utils.R"))
source(paste0(WORK_DIR,"/src/LTN/gibbs_glasso.R"))
source(paste0(WORK_DIR,"/src/LTN/gibbs_hp_cov_mu.R"))

library(mvtnorm)
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
K = 64
# true tree is left-skewed
cnt = readRDS(paste0(cachedir,'cnt_left_lambda',lambda,'K',K,'_',id,'.rds'))

# fit LTN with balanced tree
niter = 10^4
yyl = leaf2node_special(cnt, 'balanced')
if (!file.exists(paste0(cachedir,'gibbs_balanced_lambda',lambda,'K',K,'_',id,'.rds'))){
gibbs1 = gibbs_ltn_hp_cov_mu(niter=niter,YL=yyl$leftval,Y=yyl$nodeval,SEED=1,lambda=lambda)
saveRDS(gibbs1, paste0(cachedir,'gibbs_balanced_lambda',lambda,'K',K,'_',id,'.rds'))
samp_seq=(niter/2+1):niter
MU=gibbs1$MU
OMEGA=gibbs1$OMEGA
rm(gibbs1)
mu=apply(MU[,samp_seq],1,mean)
omega=Reduce('+',OMEGA[samp_seq])/length(OMEGA[samp_seq])
sigma=solve(omega)
saveRDS(list(mu_balanced=mu,sigma_balanced=sigma), paste0(cachedir,'est_para_balanced_lambda',lambda,'K',K,'_',id,'.rds'))

}


