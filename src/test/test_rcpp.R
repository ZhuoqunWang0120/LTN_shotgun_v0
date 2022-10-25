library(Rcpp)
library(RcppArmadillo)

sourceCpp('src/test/test_rcpp.cpp')
source('src/test/test_rcpp_omega.R')


tau = matrix(1:9,3,3) 
tau = t(tau) %*% tau + diag(3)
diag(tau) = 0
omega = tau + 150*diag(3)
S = t(matrix(seq(0.2,5,length.out = 9),3,3)) %*% t(matrix(seq(0.2,5,length.out = 9),3,3))+ diag(3)
print(omega)
a = update_omega_tau(S,S/10,1,omega,tau)
print(omega)
b = update_omega_tau_r(S,S/10,1,omega,tau)



