source('src/LTN/mixed_effects_shotgun.R')
source('src/LTN/mixed_effects_shotgun_uncentered_glasso.R')
source('src/LTN/utils.R')
library(ape)
library(mvtnorm)
set.seed(1)
# set parameters
N = 200
d = 5
tree = rtree(d + 1)
taxa_abundance = c(1,1,1,1,1,1)
taxa_abundance = taxa_abundance/sum(taxa_abundance)
taxa_marker_len = c(1,1,1,1,1,1)
tlr = function(x, tree){
  yyl = count2y(x,tree)
  return(qlogis(yyl$YL/yyl$Y))
}
tlr_inv = function(x, tree){
  rt=data.tree::as.Node(tree)
  K=length(tree$tip.label)
  maxLen=max(do.call(c,(lapply(ape::nodepath(tree),length))))
  #nodepath traversal of leaves: pre-order
  nodemat=do.call(rbind,lapply(ape::nodepath(tree),function(x){if(length(x)<maxLen){c(x,rep(rt$name,maxLen-length(x)))}else{x}}))
  nodemat=apply(nodemat, 1:2, function(x){paste0('n',x)})
  rt$Set(name=1:K,traversal = 'pre-order',filterFun = data.tree::isLeaf)
  rt$Do(function(x) {
    x$leftchild=names(x$children[1])
    x$rightchild=names(x$children[2])
  }, traversal = "pre-order",
  filterFun = data.tree::isNotLeaf)
  lc=stats::na.omit(rt$Get('leftchild'))
  rc=stats::na.omit(rt$Get('rightchild'))
  nam=c(paste0('n',lc),paste0('n',rc),paste0('n',rt$name))
  return(psi2p(x, nam, nodemat))
}
tlr_marker_len = tlr(taxa_marker_len, tree)
Xtest = matrix(c(rep(0, N/2), rep(1, N/2)), ncol = 1)
Xadjust = matrix(c(rep(1,N),rnorm(N * 2)), ncol = 3)
g = 4
refflabel = rep(1:4, N/4)
gam = matrix(rnorm(d*g), ncol = g)
beta1 = matrix(rep(1, d), nrow = 1)
beta2 = matrix(rep(-1:1, d), nrow = 3, byrow = F)
var_error = diag(d)
total_count = 1000
# simulate data
psi_mean = t(apply(Xtest %*% beta1 + Xadjust %*% beta2 + t(gam[,refflabel]), 1, function(x){x + tlr_marker_len}))
psi = t(apply(psi_mean, 1, function(x){rmvnorm(1,x,var_error)}))
prob = t(apply(psi, 1, function(x){tlr_inv(x,tree)}))
cnt = t(apply(prob, 1, function(x){rmultinom(1,total_count,x)}))
yyl = apply(cnt,1,function(x){count2y(x,tree)})
Y = do.call(rbind,lapply(yyl,function(x)x$Y))
YL = do.call(rbind,lapply(yyl,function(x)x$YL))
lambda_mat = matrix(rep(10, d^2),d,d)
# diag(lambda_mat) = rep(0.00001, d)
# fit model
niter = 1000
g1 = gibbs_shotgun(N = N,
                   p = d,
                   g = g,
                   YL = YL,
                   Y = Y,
                   r = 0,
                   tlr_marker_len = tlr_marker_len,
                   Xtest = Xtest,
                   Xadjust = Xadjust,
                   adjust = T,
                   grouplabel = refflabel,
                   niter = niter,
                   reff = T,
                   gprior_m = 100,
                   reffcov = 2,
                   # lambda_fixed = lambda_mat,
                   lambda_fixed = 10,
                   verbose = T)
save.image(file = 'cache/shotgun_metagenomics/confounded_original.RData')
BETA1 = data.frame(do.call(rbind,g1$BETA1))
PMAP = apply(abs(BETA1) > .Machine$double.xmin, 2, mean)
PJAP = sum(apply(abs(BETA1) > .Machine$double.xmin, 1, sum) != 0)/nrow(BETA1)
print(PMAP)
print(PJAP)
library(ggplot2)
library(reshape2)
BETA1[,'iteration'] = 1:nrow(BETA1)
ggplot(melt(BETA1, id.vars = 'iteration')) + geom_line(aes(x = iteration, y = value)) + facet_wrap(~variable)
BETA2_mean = Reduce('+',g1$BETA2)/length(g1$BETA2)
print(round(BETA2_mean,1)) # only the coef of the column of 1's seem to be problematic
INTERCEPT = do.call(rbind,lapply(g1$BETA2, function(x){x[1,]}))
GAMMA = g1$GAM
GAMMA1 = do.call(rbind,lapply(GAMMA, function(x){x[,1]}))
MEAN1 = INTERCEPT + BETA1[,-6] + GAMMA1
matplot(INTERCEPT)
matplot(BETA1[,-6])
matplot(MEAN1, type = 'l')
matplot(GAMMA1)
OMEGA = g1$OMEGA2
omega_post = Reduce('+', OMEGA)/length(OMEGA)
round(omega_post,1)
plot(unlist(lapply(OMEGA,function(x){x[1,1]})),type='l')
# test
dimnames(GAMMA1) = list(1:nrow(GAMMA1),paste0('node ',1:5))
matplot(rownames(GAMMA1), GAMMA1, type='l', col=1:5, xlab = 'iteration', ylab = 'value')
legend('bottomright', legend=colnames(GAMMA1), 
                            pch=0.5, horiz=F, col=1:5)

# why is omega[1,1] so weird when lambda_ii are small, but normal(but biased) when lambda==10?
# check whether all lambda = 10 in the new implementation yield similar result as the previous one
