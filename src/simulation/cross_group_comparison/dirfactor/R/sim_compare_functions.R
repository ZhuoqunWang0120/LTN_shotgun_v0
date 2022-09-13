#!/usr/bin/env Rscript

# source("/work/zw122/pg_tree/cross_group/dirfactor/R/MCMC.R")
# source("/work/zw122/pg_tree/cross_group/dirfactor/R/utilities.R")
library( MASS )
library( truncnorm )
library( mvtnorm )
library( MCMCpack )
library( mvtnorm )
library( tmvtnorm )
nod.free.mcmc.new = function( data, start, hyper,
                              sigma.value, sigma.prior,
                              save_path,
                              save.obj = c("sigma", "Q", "T.aug", "X", "Y.sub", "delta", "phi", "x", "er", "x.sigma" ),
                              burnin = 0.2, thin = 5, step = 1000 ){
  nv = hyper$nv
  a.er = hyper$a.er
  b.er = hyper$b.er
  #er is the variance of the pure error
  a1 = hyper$a1
  a2 = hyper$a2
  m = hyper$m
  #x.sigma variance of the regression coef
  a.x.sigma = hyper$a.x.sigma
  b.x.sigma = hyper$b.x.sigma
  #design matrix for subject
  #dim = n.sub*n.sample
  sub.design = hyper$sub.design
  sub.rep = rowSums(sub.design)

  n = ncol( data )
  n.sub = nrow( sub.design )
  p = nrow( data )
  sum.species = rowSums( data )
  sum.sample = colSums( data )
  sv.log = log( sigma.value )
  sp.log = log( sigma.prior )

  #fixed effect
  x.old = start$x
  x.sigma.old = start$x.sigma
  y.fix = start$y

  #random effect
  sigma.old = as.vector( start$sigma )
  Q.old = start$Q
  T.aug.old = as.vector( start$T.aug )
  er.old = start$er
  X.old = start$X
  Y.sub.old = start$Y.sub
  Y.old = Y.sub.old%*%sub.design

  delta.old = as.vector(start$delta)
  phi.old = start$phi

  all.cache = c()
  t.t = proc.time()
  for( iter in 2:step ){
    # browser()
    if( iter %% 1000 == 0 )
      print( iter )
    #sample sigma
    #this is a major time consuming step
    all.cache$sigma = sigma.old
    #     set.seed(1)

    sigma.old = sigma.gibbs(
      sigma.old, T.aug.old,
      Q.old, sum.species,
      sigma.value, sigma.prior,
      sv.log, sp.log
    )

    #     Q.pos = Q.old*(Q.old>0);
    #     tmp.weights = Q.pos^2%*%T.aug.old;
    #     sigma.old.t = sigma_gibbs( length( sigma.value ),
    #                              p,
    #                              sum.species,
    #                              tmp.weights,
    #                              sigma.value,
    #                              sigma.prior )


    #sample T
    all.cache$T.aug = T.aug.old
    T.aug.old = T.aug.gibbs(
      T.aug.old, sigma.old,
      Q.old, sum.sample
    )

    #sample Q
    #matrix of all parameters
    #notice the mean of the prior conditional on fixed effect and random effect is changed
    Q.para.all = cbind( c(data), c(t(X.old)%*%Y.old+t(x.old)%*%y.fix),
                        rep( sigma.old, n ), c(Q.old),
                        rep( T.aug.old, each = p ), rep( sqrt(er.old), n*p ) )
    all.cache$Q = Q.old
    #get samples

    labels = (Q.para.all[,1]==0)
    Q.0 = Q.gibbs.vec.0( Q.para.all[labels, c(2,3,5,6)] )
    Q.n0 = Q.gibbs.vec.n0( Q.para.all[!labels, c(1,2,3,5,6)],
                           Q.para.all[!labels, 4] )
    #save it into proper locations
    Q.tmp = rep( 0, nrow( Q.para.all ) )
    Q.tmp[labels] = Q.0
    Q.tmp[!labels] = Q.n0
    #convert it back to matrix
    Q.old = matrix( Q.tmp, nrow = p )
    # browser()

    #sample Y.sub
    #this is a major time consuming step
    #get tau
    tau.tmp = cumprod( delta.old )
    all.cache$Y.sub = Y.sub.old
    Y.mean.pre = X.old%*%(Q.old-t(x.old)%*%y.fix)%*%t(sub.design)/er.old

    Y.sub.old = vapply( 1:n.sub, function(x){
      Sigma.Y = solve( chol( sub.rep[x]*X.old%*%t(X.old)/er.old + diag( phi.old[x,]*tau.tmp ) ) )
      mu.Y = t(Sigma.Y)%*%Y.mean.pre[,x,drop=F]
      Sigma.Y%*%( rnorm( length(mu.Y) ) + mu.Y )
    }, rep( 0, nrow( Y.sub.old ) ) )
    Y.old = Y.sub.old%*%sub.design
    #     Y.old = matrix( 0, nrow = nrow(Y.old), ncol = ncol( Y.old ) )

    #sample er
    #prior is ga(1,10)
    all.cache$er = er.old
    er.old = 1/rgamma( 1, shape = n*p/2+a.er,
                       rate = sum((Q.old-t(X.old)%*%Y.old-t(x.old)%*%y.fix)^2)/2+b.er )

    #sample X
    Sigma.X = solve( Y.old%*%t(Y.old)/er.old + diag( m ) )
    all.cache$X = X.old
    #     X.mean = vapply( 1:p, function(x){
    #       Sigma.X%*%Y.old%*%Q.old[x,]/er.old
    #     }, rep( 0, nrow( X.old ) ) )
    X.mean = t(Sigma.X)%*%Y.old%*%t(Q.old-t(x.old)%*%y.fix)/er.old
    X.old = X.mean + t( rmvnorm( p, sigma = Sigma.X ) )

    #sample phi
    #prior is ga(nv/2,nv/2)
    #nv = 3
    all.cache$phi = phi.old
    phi.old = matrix( rgamma( n.sub*m, shape = (nv+1)/2, rate = c( Y.sub.old^2*tau.tmp + nv )/2 ),
                      nrow = n.sub, byrow = T )

    #sample delta
    #we choose a1 = 3, a2 = 4
    delta.tmp = delta.old
    all.cache$delta = delta.old
    for( k in 1:m ){
      if( k == 1 ){
        delta.prime = cumprod( delta.tmp )/delta.tmp[1]
        delta.tmp[k] = rgamma( 1, shape = n.sub*m/2 + a1,
                               rate = 1 + sum(colSums( t(Y.sub.old^2)*phi.old )*delta.prime)/2
        )
      }
      else{
        delta.prime = cumprod( delta.tmp )[-(1:(k-1))]/delta.tmp[k]
        delta.tmp[k] = rgamma( 1, shape = n.sub*(m-k+1)/2 + a2,
                               rate = 1 + sum( colSums( t(Y.sub.old^2)*phi.old )[-(1:(k-1))]*delta.prime)/2
        )
      }
    }
    delta.old = delta.tmp

    #sample the fixed effect reg coef variances
    all.cache$x.sigma = x.sigma.old
    a.x.sigma.use = a.x.sigma + p/2
    b.x.sigma.use = b.x.sigma + rowSums(x.old^2)/2
    x.sigma.old = 1/rgamma(length(b.x.sigma.use), shape = a.x.sigma.use, rate = b.x.sigma.use)

    #sample the species specific effect x
    #x is a matrix, dimension = K*I, where K is the number of covariates
    all.cache$x = x.old
    x.cov = solve( (y.fix%*%t(y.fix))/er.old + diag(1/x.sigma.old, nrow=length(x.sigma.old), ncol=length(x.sigma.old)) )
    x.mean.all = x.cov%*%y.fix%*%t(Q.old-t(X.old)%*%Y.old)/er.old
    x.old = x.mean.all + t( rmvnorm( ncol(x.mean.all), sigma = x.cov ) )

    if( iter > burnin*step & iter%%thin == 0 )
      saveRDS( all.cache[save.obj], file = paste( save_path, iter-1, sep = "/") )
  }

  #save last run
  all.cache = list( sigma = sigma.old, T.aug = T.aug.old, Q = Q.old, Y.sub = Y.sub.old,
                    X = X.old, phi = phi.old, delta = delta.old, x = x.old,
                    x.sigma = x.sigma.old, er = er.old )
  saveRDS( all.cache[save.obj], file = paste( save_path, iter, sep = "/") )

  print( "Total Time:")
  print( proc.time()-t.t )
  #   list( Q = Q, T.aug = T.aug, sigma = sigma, Y = Y, er = er, X = X, phi = phi, delta = delta )
}


