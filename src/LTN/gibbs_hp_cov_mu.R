
#' Gibbs sampler for posterior inference of LTN with sparse Gaussian graph of nodes
#' @param niter number of Gibbs iteration
#' @param Y d x N matrices of y(A)
#' @param YL d x N matrices of y(Al)
#' @param r hyperparameter for shrinkage parameter lambda. lambda ~ Ga(r,s).
#' @param s hyperparameter for shrinkage parameter lambda. lambda ~ Ga(r,s).
#' @param SEED random seed for initializing the parameters
#' @param lambda 'hp' -- the Gamma hyperprior will be used; or a real number used as a fixed lambda, and r and s will be ignored
#' @param hp_cov_mu prior covariance for mu
#' @param verbose if TRUE, print number of iterations

gibbs_glasso_hp_cov_mu=function(niter,YL,Y,r=1,s=0.01,SEED=1,lambda='hp',hp_cov_mu=NULL,verbose=F){
  #initialization
  set.seed(SEED)
  LAM=rep(0,niter)
  if (lambda=='hp'){
    lambda=stats::rgamma(1,r,s)
  }
  LAM[1]=lambda
  nsim=ncol(Y)
  kappa=NULL
  for (i in 1:nsim){
    kappa=cbind(kappa,YL[,i]-Y[,i]/2)
  }
  p=nrow(Y)
  K=p+1
  PHI=diag(rep(1,p))
  OMEGA=list()
  set.seed(SEED)
  OMEGA[[1]]=stats::rWishart(1,p+2,PHI)[,,1]
  omega=OMEGA[[1]]
  TAU=list()
  tau=matrix(0,nrow=p,ncol=p)
  set.seed(SEED)
  for (l in 2:p){
    for (k in 1:(l-1)){
      mu_prime=sqrt(lambda^2/omega[k,l]^2)
      ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
      tau[k,l]=1/ukl
      tau[l,k]=1/ukl
    }
  }
  TAU[[1]]=tau
  PSI=list()
  set.seed(SEED)
  psi=matrix(0,ncol = ncol(Y),nrow=nrow(Y))
  for (i in 1:nrow(Y)){
    for (j in 1:ncol(Y)){
      if (Y[i,j]==0){
        psi[i,j]=0
      }
      else{
        tmp=YL[i,j]/Y[i,j]
        if (tmp==0){
          psi[i,j]=-5
        }
        else{
          if (tmp==1){
            psi[i,j]=5
          }
          else{
            psi[i,j]=stats::qlogis(tmp)
          }
        }
      }
    }
  }
  PSI[[1]]=psi
  if (is.null(hp_cov_mu)){
    Lam = diag(rep(5,p))
  }
  else{
    Lam = hp_cov_mu
  }
  Z=list()
  z=matrix(0,ncol=nsim,nrow=p)
  set.seed(SEED)
  for (i in 1:nsim){
    for (j in 1:p){
      if (Y[j,i]!=0){
        z[j,i]=BayesLogit::rpg(1,Y[j,i]+0.0001,0)
        while (is.na(z[j,i])){
          z[j,i]=BayesLogit::rpg(1,Y[j,i]+0.0001,0)
          warning('NA in PG variable')
        }
      }
    }
  }
  Z[[1]]=z
  set.seed(SEED)
  MU=matrix(0,ncol = niter,nrow=p)
  #gibbs sampling
  for (it in 2:niter){
    if (verbose){print(it)}
    # print(it)
    lambda=LAM[it-1]
    omega=OMEGA[[it-1]]
    psi=PSI[[it-1]]
    tau=TAU[[it-1]]
    mu=MU[,it-1]
    #update z
    z=matrix(0,ncol=nsim,nrow=p)
##
z_na = list()
##
    for (i in 1:nsim){
      for (j in 1:p){
        if (Y[j,i]!=0){ #only update nodes with non-zero y_i(A)
          #z[j,i]=BayesLogit::rpg(1,as.numeric(Y[j,i]),psi[j,i])
		z[j,i] = rpg0(Y[j,i],psi[j,i])
		##
if (is.na(z[j,i])){
z_na = c(z_na,list(j=j,i=i,z=z[j,i],para1=as.numeric(Y[j,i]),para2=psi[j,i]))
}
##

        }
      }
    }
    Z[[it]]=z
    #update psi
    psi=matrix(0,ncol = nsim,nrow=p)
    for (i in 1:nsim){
      # omega errors
      t1 <- try(cov_mat<-chol2inv(chol(omega+diag(z[,i]))))
  
##
    if ("try-error" %in% class(t1)){
        # break
	warning('omega+z not invertible')
        return(list(iteration=it,OMEGA=OMEGA, omega_it=omega,z_it=z,Z=Z,psi_it=psi,MU=MU,z_na=z_na))
      }
##


    mean_vec=cov_mat%*%(omega%*%mu+kappa[,i])
      psi[,i]=mvtnorm::rmvnorm(1,mean_vec,cov_mat)
    }
    PSI[[it]]=psi
    #update mu
    cov_mat=chol2inv(chol(chol2inv(chol(Lam))+nsim*omega))
    mean_vec=cov_mat%*%omega%*%apply(psi,1,sum)
    MU[,it]=mvtnorm::rmvnorm(1,mean_vec,cov_mat)
    mu=MU[,it]
    #update omega and tau together (blocked gibbs)
    psi_diff=psi-matrix(rep(mu,nsim),ncol=nsim,byrow=F)
    S=psi_diff%*%t(psi_diff)
    #step1
    for (k in 1:p){
      #(b)
      s22=S[k,k]
      gam=stats::rgamma(1,nsim/2+1,(s22+lambda)/2)
      Dtau=diag(tau[-k,k])
      omega11=omega[-k,-k]
      C=chol2inv(chol((s22+lambda)*chol2inv(chol(omega11))+chol2inv(chol(Dtau))))
      s21=S[-k,k]
      beta_mean=-C%*%s21
      be=mvtnorm::rmvnorm(1,beta_mean,C)
      #(c)
      omega[-k,k]=be
      omega[k,-k]=t(be)
      omega[k,k]=gam+be%*%chol2inv(chol(omega[-k,-k]))%*%t(be)
    }
    # step2
    for (l in 2:p){
      for (k in 1:(l-1)){
        mu_prime=sqrt(lambda^2/omega[k,l]^2)
        ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
        tau[k,l]=1/ukl
        tau[l,k]=1/ukl
      }
    }
    OMEGA[[it]]=omega
    TAU[[it]]=tau
    if (lambda=='hp'){
      # update lambda
      lambda=stats::rgamma(1,r+p*(p+1)/2,s+sum(abs(omega))/2)
    }
    LAM[it]=lambda
  }
  return(list(OMEGA=OMEGA,LAM=LAM,MU=MU))
}

#' Gibbs sampler for posterior inference of LTN with sparse Gaussian graph of nodes
#' @param data phyloseq object containing OTU table and phylogenetic tree
#' @param Y (numeric) n*d matrix of y_i(A). Y and YL will be used if the phyloseq object is not provided
#' @param YL (numeric) n*d matrix of y_i(A_l)
#' @param niter number of Gibbs iterations
#' @param SEED random seed for initializing the parameters
#' @param lambda shrinkage parameter in graphical lasso prior
#' @param hp_cov_mu prior covariance of mu
#' @param verbose if TRUE, print number of iterations
#' @export
gibbs_ltn_hp_cov_mu=function(data=NULL,Y=NULL,YL=NULL,lambda=10,niter,SEED=1,hp_cov_mu=NULL,verbose=F){
  if (!is.null(data)){
    # compute y(A) and y(Al)
    tree=phyloseq::phy_tree(data)
    cnt=phyloseq::otu_table(data)
    if (phyloseq::taxa_are_rows(data)){
      cnt=t(cnt)
    }
    cnt=cnt[,tree$tip.label]
    yyl=apply(cnt,1, function(x) {count2y(x,tree)})
    Y=do.call(rbind,lapply(yyl, function(x){x[['Y']]}))
    YL=do.call(rbind,lapply(yyl, function(x){x[['YL']]}))
  }
  return(gibbs_glasso_hp_cov_mu(niter=niter,YL=t(YL),Y=t(Y),r=1,s=0.01,SEED=SEED,lambda=lambda,hp_cov_mu=hp_cov_mu,verbose=verbose))
}

rpg0 = function(h, z) {
  h = as.numeric(h)
  z = as.numeric(z)
  if (h == 0) {
    x = 0
  }
  else{
    x = BayesLogit::rpg(1, h, z)
  }
  while (is.na(x)) {
    if (h == 0) {
      x = 0
    }
    else{
      x = BayesLogit::rpg.gamma(1, h, z)
    }
  }
  return(x)
}



