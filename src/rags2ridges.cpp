// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// To avoid "::" usage uncomment below:
//using namespace Rcpp;
//using namespace RcppArmadillo;
//using namespace arma;


/* -----------------------------------------------------------------------------

   AUXILIARY TOOLS

----------------------------------------------------------------------------- */



// [[Rcpp::export]]
arma::mat armaPooledS(const Rcpp::List & SList,  // List of covariance matrices
                      const Rcpp::NumericVector ns,
                      const int mle = 0) {
  /* ---------------------------------------------------------------------------
   Function to compute the pooled covariance estimate. Returns a numeric matrix
   giving the pooled covariance matrix.
   - SList > A list of covariance matrices.
   - nu    > A numeric vector giving the number of samples corresponding
             to each scatter matrix.
   - mle   > A integer of length one equalling 0 or 1. If 0, the bias
             corrected ML is used. If 1 the ML estimate is used.
  --------------------------------------------------------------------------- */

  const int K = SList.size();
  const int imle = 1 - mle;
  const double fac = 1.0f/(sum(ns) - K*imle);
  arma::mat S0 = SList[0];
  S0 = (ns[0] - imle)*S0;
  for (int i = 1; i < K; ++i) {
    arma::mat Si = SList[i];
    S0 += (ns[i] - imle)*Si;
  }
  return fac*S0;
}



/* -----------------------------------------------------------------------------

   GENERAL RIDGE AND FUSED RIDGE ESTIMATION TOOLS

----------------------------------------------------------------------------- */



arma::mat armaEigShrink(const arma::vec eigvals,
                        const double lambda,
                        const double alpha) {
  /* ---------------------------------------------------------------------------
   Function that shrinks the eigenvalues
   Shrinkage is that of the rotation equivariant alternative ridge estimator
   - eigvals is vector on eigenvalues
   - lambda is a double giving the penalty
   - alpha is the unique diagonal value in the rotaional equivariant target
  --------------------------------------------------------------------------- */


  arma::vec seigvals = eigvals - lambda*alpha;
  return sqrt(lambda + 0.25f*pow(seigvals, 2.0f)) + 0.5f*seigvals;
}



arma::mat armaRidgeSAnyTarget_OLD(const arma::mat & S,
                                 const arma::mat & target,
                                 const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general targets.
   - S is the sample covariance matrix
   - target is the target matrix with the same size as S
   - lambda is the ridge penalty
  --------------------------------------------------------------------------- */

  const int n = S.n_cols;
  const arma::mat E = symmatl(S) - lambda * symmatl(target);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, 0.25f*E*E + lambda*arma::eye(n, n), "dc");

  return inv_sympd(0.5f*E + eigvec*diagmat(sqrt(eigval))*eigvec.t());
}



arma::mat armaRidgeSAnyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general targets.
   - S is the sample covariance matrix
   - target is the target matrix with the same size as S
   - lambda is the ridge penalty
  --------------------------------------------------------------------------- */

  const int n = S.n_cols;
  const double inv_lambda = 1.0f/lambda;
  const arma::mat E = symmatl(S) - lambda * symmatl(target);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, 0.25f*E*E + lambda*arma::eye(n, n), "dc");

  eigval = inv_lambda*sqrt(eigval);
  return (-0.5f*inv_lambda)*E + eigvec*diagmat(eigval)*eigvec.t();
}



arma::mat armaRidgeSRotationInvariantTarget_OLD(const arma::mat & S,
                                                const double alpha,
                                                const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for rotational equivariate target.
   - S is the sample covariance matrix
   - alpha is a scaling of the identity matrix
   - lambda is the ridge penalty
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, symmatl(S), "dc");  // Eigen decomposition

  eigvals = armaEigShrink(eigvals, lambda, alpha);

  return inv_sympd(eigvecs*diagmat(eigvals)*eigvecs.t());
}



arma::mat armaRidgeSRotationInvariantTarget(const arma::mat & S,
                                            const double alpha,
                                            const double lambda) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in the rotational equivariate case
   - S is the sample covariance matrix
   - alpha is a scaling of the identity matrix
   - lambda is the ridge penalty
  --------------------------------------------------------------------------- */

  const double inv_lambda = 1.0f/lambda;
  arma::mat E = symmatl(S);
  E.diag() -= alpha*lambda;

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, symmatl(S), "dc");  // Eigen decomposition

  eigvals = armaEigShrink(eigvals, lambda, alpha);

  return inv_lambda*(symmatl(eigvecs*diagmat(eigvals)*eigvecs.t())-symmatl(E));
}



// [[Rcpp::export]]
arma::mat armaRidgeS(const arma::mat & S,
                     const arma::mat & target,
                     const double lambda) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in C++.
   - S      > the sample covariance matrix (a numeric matrix on the R side)
   - target > target matrix (a numeric matrix on the R side, same size as S)
   - lambda > the penalty (a numeric of length one on the R side)
  --------------------------------------------------------------------------- */

  const int n = S.n_rows;
  const double alpha = target(0, 0);
  const arma::mat alphaI = alpha*arma::eye<arma::mat>(n, n);
  if (arma::all(arma::all(target == alphaI))) {
    return armaRidgeSRotationInvariantTarget(S, alpha, lambda);
  } else {
    return armaRidgeSAnyTarget(S, target, lambda);
  }
}



// [[Rcpp::export]]
arma::mat armaFusedUpdate(int k0,
                          const Rcpp::List & PList,
                          const Rcpp::List & SList,
                          const Rcpp::List & TList,
                          const arma::vec ns,
                          const double lambda,
                          arma::mat lambdaFmat) {
  /* ---------------------------------------------------------------------------
   "Update" the covariance matrices and use the regular ridge estimate.
   - k0          > An integer giving the class estimate to be updated.
   - PList       > A list of length K of matrices giving the current precision
                   estimates.
   - SList       > A list of length K of sample correlation matrices the same
                   size as those of PList.
   - TList       > A list of length K of target matrices the same size
                   as those of PList
   - ns          > A vector of length K giving the sample sizes.
   - lambda      > The ridge penalty (a postive number).
   - lambdaFmat  > A K by K symmetric adjacency matrix giving the fused
                   penalty graph with non-negative entries where
                   lambdaFmat[k1, k2] determine the (rate of) shrinkage
                   between estimates in classes corresponding to SList[k1]
                   and SList[k1].
  --------------------------------------------------------------------------- */

  k0 = k0 - 1;  // Shift index to C++ convention
  const int n = ns.n_elem;
  const int K = SList.size();
  lambdaFmat(k0, k0) = 0;
  const double a = (sum(lambdaFmat.row(k0)) + lambda)/(ns[k0]);

  arma::mat S0 = Rcpp::as<arma::mat>(Rcpp::wrap(SList[k0]));
  arma::mat T0 = Rcpp::as<arma::mat>(Rcpp::wrap(TList[k0]));
  arma::mat O(n, n);
  arma::mat T(n, n);
  for (int k = 0; k < K; ++k) {
     if (k == k0) {
       continue;
     }
     O = Rcpp::as<arma::mat>(Rcpp::wrap(PList[k]));
     T = Rcpp::as<arma::mat>(Rcpp::wrap(TList[k]));
     S0 -= (lambdaFmat(k, k0)/ns(k0))*(O - T);
  }
  return armaRidgeS(S0, T0, a);
}



/* -----------------------------------------------------------------------------

   GENERAL SIMULATION TOOLS

----------------------------------------------------------------------------- */



// [[Rcpp::export]]
arma::mat rmvnormal(const int n, arma::rowvec mu, arma::mat sigma) {
  /* ---------------------------------------------------------------------------
   Simulate from multivariate normal distribution with mean mu and covariance
   sigma. Returns a matrix of size n by length(mu) of observations.
   (Code taken from package GMCM)
   - n     > An integer giving the number of samples to simulate.
   - mu    > A vector giving the population mean.
   - sigma > A matrix giving the population covariance matrix.
  --------------------------------------------------------------------------- */

  Rcpp::RNGScope();  // Allows for using set.seed(...) on the R side
  const int d = mu.size();

  // Create matrix of standard normal random values
  arma::mat ans(n, d, arma::fill::none);
  for (int j = 0; j < d; ++j) {  // Fill ans with random values
    ans.col(j) = Rcpp::as<arma::colvec>(Rcpp::rnorm(n));
  }

  // Do the Cholesky decomposition
  const arma::mat csigma = arma::chol(sigma);

  // Do the transformation
  ans = ans * csigma;
  ans.each_row() += mu; // Add mu to each row in transformed ans

  return ans;
}



arma::mat armaRWishartSingle(const arma::mat L, const double nu, const int p) {
  // Wishart simulation capabilities: (taken from package correlateR)
  // Simualte a single wishart distributed matrix
  arma::mat A(p, p, arma::fill::randn);
  A = trimatl(A);
  arma::vec chisq(p);
  for (int i = 0; i < p; ++i) {
    chisq[i] = sqrt(Rcpp::as<double>(Rcpp::rchisq(1, nu - i)));
  }
  A.diag() = chisq;
  arma::mat LA = L*A;
  return LA.t()*LA;
}

// [[Rcpp::export]]
arma::cube armaRWishart(const int n,
                        const arma::mat & sigma,
                        const double nu) {
  // Simulate n wishart distributed matrices
  const int p = sigma.n_cols;
  const arma::mat L = arma::chol(sigma);
  arma::cube ans(p, p, n);
  for (int k = 0; k < n; ++k) {
    ans.slice(k) = armaRWishartSingle(L, nu, p);
  }
  return ans;
}

// [[Rcpp::export]]
arma::cube armaRInvWishart(const int n,
                           const arma::mat & psi,
                           const double nu) {
  // Simulate n inverse wishart distributed matrices
  const int p = psi.n_cols ;
  const arma::mat L = arma::chol(inv(psi));
  arma::cube ans(p, p, n);
  for (int k = 0; k < n; ++k) {
    ans.slice(k) = inv(armaRWishartSingle(L, nu, p));
  }
  return ans;
}



/*** R
# Testing
library("rags2ridges")
library("microbenchmark")
S  <- createS(n = 10, p = 10)
target <- default.target(S, type = "DEPV")
lambda <- 2



# General target
ridgeS1 <- function(S, target, lambda) {
  E <- (S - lambda * target)
  return(solve(E/2 + expm::sqrtm((E %*% E)/4 + lambda * diag(nrow(S)))))
}

microbenchmark(A1 <- ridgeS1(S, target, lambda),
               B1 <- armaRidgeSAnyTarget(S, target, lambda),
               C1 <- armaRidgeS(S, target, lambda),
               D1 <- ridgeS(S, lambda, target = target),
               times = 1)
stopifnot(all.equal(unname(A1), B1))
stopifnot(all.equal(unname(A1), C1))

# NULL TARGET
ridgeS2 <- function(S, lambda) {
  Spectral  <- eigen(S, symmetric = TRUE)
  Eigshrink <- rags2ridges:::.eigShrink(Spectral$values, lambda)
  return(solve(Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)))
}
target <- default.target(S, "Null")
microbenchmark(A2 <- ridgeS2(S, lambda),
               B2 <- armaRidgeSRotationInvariantTarget(S, 0.0, lambda),
               C2 <- armaRidgeS(S, target, lambda))
stopifnot(all.equal(unname(A2), unname(B2)))
stopifnot(all.equal(unname(A2), unname(C2)))

# Equal diagnonal target
ridgeS3 <- function(S, target, lambda) {
  varPhi    <- unique(diag(target))
  Spectral  <- eigen(S, symmetric = TRUE)
  Eigshrink <- rags2ridges:::.eigShrink(Spectral$values, lambda, const = varPhi)
  return(solve(Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)))
}

target <- default.target(S) # Equal diagnoal
microbenchmark(A3 <- ridgeS3(S, target, lambda),
               B3 <- armaRidgeSRotationInvariantTarget(S, target[1,1], lambda),
               C3 <- armaRidgeS(S, target, lambda))
stopifnot(all.equal(unname(A3), unname(B3)))
stopifnot(all.equal(unname(A3), unname(C3)))



library("rags2ridges")
library("microbenchmark")
S  <- createS(n = 10, p = 500)
target <- default.target(S, type = "DEPV")#default.target(S)#
lambda <- 2
system.time(B <- ridgeSArma(S, target = target, lambda = lambda))
microbenchmark(A <- ridgeS(S, target = target, lambda = lambda),
               B <- armaRidgeS(S, target = target, lambda = lambda),
               times = 1)
stopifnot(all.equal(A, B))

S  <- createS(n = 10, p = 5000)
target <- default.target(S, type = "DEPV")
lambda <- 2
system.time(B <- ridgeSArma(S, target = target, lambda = lambda))
# Takes about 2 mins -- about 10-20 hours for the ridgeS


#
# Test armaFusedUpdate
#
ns <- c(100, 100)
PList <- createS(n = c(100, 100), p = 10)
SList <- createS(n = c(100, 100), p = 10)
TList <- default.target.fused(SList, ns)
k0 <- 1

res <- microbenchmark(
  A = armaFusedUpdate2(k0, PList, SList, TList, ns, 1, matrix(1, 2, 2)),
  B = armaFusedUpdate(k0, PList, SList, TList, ns, 1, matrix(1, 2, 2)),
  C = rags2ridges:::.fusedUpdate(k0, PList, SList, TList, ns, 1, matrix(1, 2, 2)),
  times = 10000)
boxplot(res)

all.equal(A, C)
all.equal(A2, A)
all.equal(A, B)
*/

