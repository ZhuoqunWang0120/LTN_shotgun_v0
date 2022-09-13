#!/usr/bin/env Rscript

argv=commandArgs(TRUE)
id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
i = ceiling(id / 2)
n_or_a = 2 * i - id
WORK_DIR=argv[1]
r=as.numeric(argv[2])
niter=as.numeric(argv[3])
workdir=paste0(WORK_DIR,'/src/simulation/cross_group_comparison/dirfactor/')
source(paste0(workdir,"/R/MCMC.R"))
source(paste0(workdir,"/R/utilities.R"))
source(paste0(workdir,'/R/sim_compare_functions.R'))
resdir=paste0(WORK_DIR,'/cache/cross_group_comparison/multi_otu/dirfactor/')
system(paste0('mkdir -p ',resdir,'/gibbs'))
system(paste0('mkdir -p ',resdir,'/teststat'))
data1dir=paste0(WORK_DIR,"/cache/ps_sim_100.RData")
data2dir=paste0(WORK_DIR,"/cache/cross_group_comparison/multi_otu/")
if (!file.exists(paste0(resdir,'teststat/','i',i,'H',n_or_a,'r',r,'iter',niter,'_','v_post_norm.RData'))){
	input_data1=readRDS(data1dir)
	input_data2=readRDS(paste0(data2dir,'sim',i,'H',n_or_a,'.RData'))
	# data
	cnt=input_data2$cnt
	data=t(cnt)
	n = ncol(data)
	p = nrow(data)
	# blocking factors (random effects)
	grouplabel=data.frame(id=as.character(input_data1$grouplabel))
	sub.design = t(model.matrix( ~.-1, grouplabel ))
	n.sub = nrow(sub.design)
	# hyper
	# m is number of factors
	# all hyperparameters other than m=r are default
	hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = r,
		     sub.design = sub.design, alpha = 10, beta = 0, a.x.sigma=5, b.x.sigma=5 )
	sigma.value = seq(0.001,0.999,0.001)
	tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
	sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
	sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
	# fixed effects
	Xtest=input_data2$Xtest
	Xadjust=input_data1$Xadjust
	y.fix.inter=t(cbind(Xtest,Xadjust))
	# permutation
	# y.fix.inter[1,] = matrix(sample(y.fix.inter[1,]),nrow = 1) #shuffle
	# start
	start.sub = list( er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
			                   sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ),
					                     T.aug = rgamma( n, 10, 1 ),
					                     Q = matrix( 0.5, nrow = p, ncol = n ),
							                       X = matrix( rnorm( hyper$m*p ), nrow = hyper$m ),
							                       Y.sub = matrix( 0, nrow = hyper$m, ncol = n.sub ),
									                         delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ), rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
									                         phi = matrix( rgamma( hyper$m*n.sub, shape = 3/2, rate = 3/2 ), nrow = n.sub ),
												                   y = y.fix.inter,
												                   x = matrix( rnorm( p*nrow(y.fix.inter) ), nrow=nrow(y.fix.inter) ),
														                     x.sigma = rgamma( nrow(y.fix.inter), 1, 1 ) )
	# mkdir
	# mkdir
	if(!dir.exists(paste0(resdir,'/gibbs/i',i,'H',n_or_a,'r',r,'iter',niter))){
		  dir.create(paste0(resdir,'/gibbs/i',i,'H',n_or_a,'r',r,'iter',niter), recursive = T)
	}
	# mcmc
	nod.free.mcmc.new( data, start.sub, hyper, sigma.value, sigma.prior, paste0(resdir,'/gibbs/i',i,'H',n_or_a,'r',r,'iter',niter), burnin = 0.5, step = niter, thin = 10 )

	# calculate teststat
	sampdir=paste0(resdir,'/gibbs/i',i,'H',n_or_a,'r',r,'iter',niter,'/')
	V=sapply(list.files(sampdir),function(filename){readRDS(paste0(sampdir,filename))$x[1,]})
	v_post=apply(V,1,mean)
	v_post_norm=sqrt(sum(v_post^2))
	saveRDS(v_post,paste0(resdir,'teststat/','i',i,'H',n_or_a,'r',r,'iter',niter,'_','v_post.RData'))
	saveRDS(v_post_norm,paste0(resdir,'teststat/','i',i,'H',n_or_a,'r',r,'iter',niter,'_','v_post_norm.RData'))
}


