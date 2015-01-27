// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// These are not needed:
//using namespace Rcpp;
//using namespace RcppArmadillo;
//// [[Rcpp::depends(RcppArmadillo)]]


arma::mat armaEigShrink(const arma::vec eigvals,
                        const double lambda,
                        const double constant) {
  // Function that shrinks the eigenvalues in an eigenvector
  // Shrinkage is that of the rotation equivariant alternative ridge estimator
  arma::vec seigvals = eigvals - lambda*constant;
  return sqrt(lambda + 0.25f*pow(seigvals, 2.0f)) + 0.5f*seigvals;
}


// [[Rcpp::export]]
arma::mat armaRidgeSAnyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda) {
  // S is the sample covariance matrix
  // target is the target matrix with the same size as S
  // lambda is the ridge penalty
  const int n = S.n_cols;

  const arma::mat E = symmatl(S) - lambda * symmatl(target);
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, 0.25f*E*E + lambda*arma::eye(n, n), "dc");
  return inv_sympd(0.5f*E + eigvec*diagmat(sqrt(eigval))*eigvec.t());
}

// [[Rcpp::export]]
arma::mat armaRidgeSZeroTarget(const arma::mat & S,
                               const double lambda) {
  // S is the sample covariance matrix
  // lambda is the ridge penalty
  // Here, the target is the zero matrix (and hence not explicitly needed)
  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, S, "dc");  // Eigen decomposition

  eigvals = armaEigShrink(eigvals, lambda, 0.0f);

  return inv_sympd(eigvecs*diagmat(eigvals)*eigvecs.t());
}

// [[Rcpp::export]]
arma::mat armaRidgeSEqualDiagTarget(const arma::mat & S,
                                    const arma::mat & target,
                                    const double lambda) {
  // S is the sample covariance matrix
  // target is the zero target matrix
  // lambda is the ridge penalty
  arma::vec var_phi = unique(target.diag());
  if (var_phi.n_elem != 1) {
    throw std::runtime_error("Number of unique elements in diagnoal is not 1");
  }

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, S, "dc");  // Eigen decomposition

  eigvals = armaEigShrink(eigvals, lambda, var_phi(0));

  return inv_sympd(eigvecs*diagmat(eigvals)*eigvecs.t());
}


/*** R
# Testing
library("rags2ridges")
library("microbenchmark")
S  <- createS(n = 10, p = 50)
target <- default.target(S, type = "DEPV")
lambda <- 2


# General target
ridgeS1 <- function(S, target, lambda) {
  E <- (S - lambda * target)
  return(solve(E/2 + expm::sqrtm((E %*% E)/4 + lambda * diag(nrow(S)))))
}

microbenchmark(A1 <- ridgeS1(S, target, lambda),
               B1 <- armaRidgeSAnyTarget(S, target, lambda))
stopifnot(all.equal(unname(A1), B1))


# NULL TARGET
ridgeS2 <- function(S, lambda) {
  Spectral  <- eigen(S, symmetric = TRUE)
  Eigshrink <- rags2ridges:::.eigShrink(Spectral$values, lambda)
  return(solve(Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)))
}
microbenchmark(A2 <- ridgeS2(S, lambda),
               B2 <- armaRidgeSZeroTarget(S, lambda))
stopifnot(all.equal(unname(A2), unname(B2)))


# Equal diagnoal target
ridgeS3 <- function(S, target, lambda) {
  varPhi    <- unique(diag(target))
  Spectral  <- eigen(S, symmetric = TRUE)
  Eigshrink <- rags2ridges:::.eigShrink(Spectral$values, lambda, const = varPhi)
  return(solve(Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)))
}

target <- default.target(S) # Equal diagnoal
microbenchmark(A3 <- ridgeS3(S, target, lambda),
               B3 <- armaRidgeSEqualDiagTarget(S, target, lambda))
stopifnot(all.equal(unname(A3), unname(B3)))


library("rags2ridges")
library("microbenchmark")
S  <- createS(n = 10, p = 500)
target <- default.target(S, type = "DEPV")#default.target(S)#
lambda <- 2
system.time(B <- ridgeSArma(S, target = target, lambda = lambda))
microbenchmark(A <- ridgeS(S, target = target, lambda = lambda),
               B <- ridgeSArma(S, target = target, lambda = lambda),
               times = 1)
stopifnot(all.equal(A, B))

S  <- createS(n = 10, p = 5000)
target <- default.target(S, type = "DEPV")
lambda <- 2
system.time(B <- ridgeSArma(S, target = target, lambda = lambda)
            # Takes about 2 mins -- about 10-20 hours for the ridgeS
*/

