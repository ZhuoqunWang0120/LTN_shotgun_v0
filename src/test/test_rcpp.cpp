#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M) {
    return arma::eig_sym(M);
}

// [[Rcpp::export]]
Rcpp::List update_omega_tau(arma::mat S, arma::mat lambda, double nsim, arma::mat omega, arma::mat tau) {
    Rcpp::Environment statmod = Rcpp::Environment::namespace_env("statmod");
    Rcpp::Function rinvgaussC = statmod["rinvgauss"];
  // lambda is a matrix
  int p = S.n_rows ;
  // arma::mat omega_new(p, p) ; 
  // omega_new.fill(0.0) ; 
  arma::mat Dtau(p - 1, p - 1) ;
  Dtau.fill(0.0) ;
  arma::mat omega11(p - 1, p - 1) ;
  omega11.fill(0.0) ;
  arma::mat Cj(p - 1, p - 1) ;
  Cj.fill(0.0) ;
  arma::mat s21(p - 1,1) ;
  s21.fill(0.0) ;
  for (int k = 0; k < p; k++){
    double s22 ;
    s22 = S(k, k) ;  
    // Rcout << "The value of s22(S[k,k]) : " << s22 << "\n";
    //    gam = stats::rgamma(1, nsim / 2 + 1, (s22 + lambda[k,k]) / 2) shape-rate
    //double gam = Rcpp::rgamma(int n = 1, double shape = nsim / 2 + 1, double rate = (s22 + lambda(k - 1, k - 1)) / 2 ) ; 
    double gam = R::rgamma(nsim / 2 + 1,  (s22 + lambda(k, k)) / 2 ) ; 
    // Dtau = diag(tau[-k, k])
    int rowIndex = 0;
    for (int ii = 0; ii < p; ii++){
        if (!(ii==k)){
            Dtau(rowIndex, rowIndex) = tau(ii, k) ;
            rowIndex = rowIndex + 1 ;
        }
    }
    //Dtau.print("Dtau") ;
    // omega11 = omega[-k, -k] 
    omega11.fill(0.0) ;
    int colIndex = 0;
    rowIndex = 0;
    for (int ii = 0; ii < p; ii++){ 
        // ii is row
        if (ii!=k){
            // subset (ii, -k)
            colIndex = 0;
            for (int jj = 0; jj < p; jj++){
                if (jj!=k){
                    omega11(rowIndex,colIndex) = omega(ii,jj) ;
                    colIndex = colIndex + 1;
                }
            }
            rowIndex = rowIndex + 1 ;
        }
    }
    //omega11.print("omega11");
    // C = chol2inv(chol((s22 + lambda[k,k]) * chol2inv(chol(omega11)) + chol2inv(chol(Dtau))))
    // Cj = inv((s22 + lambda(k,k)) * inv(omega11) + inv(Dtau));
    Cj = ((s22 + lambda(k,k)) * omega11.i() + Dtau.i()).i(); 
    Cj.print("Cj") ;
    // s21 = S[-k, k]
    rowIndex = 0;
    for (int ii = 0; ii < p; ii++){
        if (!(ii==k)){
            s21(rowIndex, 0) = S(ii, k) ;
            rowIndex = rowIndex + 1 ;
        }
    }
    //s21.print("s21") ;
    // beta_mean = -C %*% s21
    arma::mat beta_mean = -Cj * s21;
    //beta_mean.print("beta_mean") ;
    // be = mvtnorm::rmvnorm(1, beta_mean, C)
    arma::mat be = arma::mvnrnd( beta_mean, Cj, 1 ); 
    //be.print("be") ;
    // (c)
    // omega[-k, k] = be
    rowIndex = 0;
    for (int ii = 0; ii < p; ii++){
        if (!(ii==k)){
            omega(ii, k) = be(rowIndex,0) ;
            rowIndex = rowIndex + 1 ;
        }
    }
    // omega[k, -k] = t(be)
    colIndex = 0;
    for (int ii = 0; ii < p; ii++){
        if (!(ii==k)){
            omega(k, ii) = be(colIndex,0) ;
            colIndex = colIndex + 1 ;
        }
    }
    omega.print("omega") ;
    // omega[k, k] = gam + be %*% chol2inv(chol(omega[-k, -k])) %*% t(be)
    // omega[-k,-k]
    colIndex = 0;
    rowIndex = 0;
    omega11.fill(0.0);
    for (int ii = 0; ii < p; ii++){ 
        // ii is row
        if (ii!=k){
            // subset (ii, -k)
            colIndex = 0;
            for (int jj = 0; jj < p; jj++){
                if (jj!=k){
                    omega11(rowIndex,colIndex) = omega(ii,jj) ;
                    colIndex = colIndex + 1;
                }
            }
            rowIndex = rowIndex + 1 ;
        }
    }
    arma::mat qprod = be.t() * omega11.i() * be;
    qprod.print("qprod") ;
    omega(k,k) = gam + qprod(0,0);
  }
  // step 2
  for (int l = 1; l < p; l++){
    for (int k = 0; k < l; k++){
        // mu_prime = sqrt(lambda[k,l] ^ 2 / omega[k, l] ^ 2)
        double mu_prime = sqrt(pow(lambda(k,l), 2)/pow(omega(k, l), 2));
        // ukl = statmod::rinvgauss(1, mu_prime, lambda[k,l] ^ 2)
        double lambda_prime = pow(lambda(k,l),2);
        Rcpp::NumericVector ukl_vec = rinvgaussC(1, mu_prime, lambda_prime); 
        double ukl_inv = 1/ukl_vec(0);
        Rprintf("the value of 1/ukl: %f \n", ukl_inv);
        tau(k, l) = ukl_inv;
        tau(l, k) = ukl_inv;
    }
  }
  Rcpp::List ret;
  ret["omega"] = omega;
  ret["tau"] = tau;
  return(ret) ;
}


// 