// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// These are not needed:
//using namespace Rcpp;
//using namespace RcppArmadillo;
//// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat RcppArmadilloRidgeS(arma::mat & S, arma::mat & target, double lambda) {
  // S is the sample covariance matrix
  // target is the target matrix with the same size as S
  // lambda is the ridge penalty
  int n = S.n_cols;

  arma::mat E = symmatl(S) - lambda * symmatl(target);
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, 0.25f*E*E + lambda*arma::eye(n, n), "dc");
  return inv_sympd(0.5f*E + eigvec*diagmat(sqrt(eigval))*eigvec.t());
}


/*** R
# Testing
library("rags2ridges")
S  <- createS(n = 10, p = 20)
target <- default.target(S, type = "DEPV")
lambda <- 2

FRest <- function(S, target, lambda) {
  E <- (S - lambda * target)
  return(solve(E/2 + expm::sqrtm((E %*% E)/4 + lambda * diag(nrow(S)))))
}

# Better benchmark
microbenchmark::microbenchmark(FRest(S, target, lambda),
                               RcppArmadilloRidgeS(S, target, lambda))

*/

