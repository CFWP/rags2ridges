// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

// To avoid "::" usage uncomment below:
//using namespace Rcpp;
//using namespace RcppArmadillo;
//using namespace arma;


////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

   AUXILIARY TOOLS

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export(.armaPooledS)]]
arma::mat armaPooledS(const Rcpp::List & Slist,  // List of covariance matrices
                      const Rcpp::NumericVector ns,
                      const int mle = 0) {
  /* ---------------------------------------------------------------------------
   Function to compute the pooled covariance estimate. Returns a numeric matrix
   giving the pooled covariance matrix.
   - Slist > A list of covariance matrices.
   - nu    > A numeric vector giving the number of samples corresponding
             to each scatter matrix.
   - mle   > A logical (or integer) of length one equalling FALSE (0) or
             TRUE (1). If FALSE, the bias corrected ML estimate is used.
             If TRUE the ML estimate is used.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  const int imle = 1 - mle;
  const double rdenum = 1.0f/(sum(ns) - G * imle);
  arma::mat S0 = Slist[0];       // First element of Slist
  S0 = (ns[0] - imle)*S0;        // Multiply by the class sample size (- 1)
  for (int i = 1; i < G; ++i) {  // Loop through remaning elements. Note i = 1.
    arma::mat Si = Slist[i];
    S0 += (ns[i] - imle)*Si;
  }
  return rdenum*S0;
}



////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

  REGULAR (NON-FUSED) RIDGE ESTIMATOR

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////


inline arma::mat armaEigShrink(const arma::vec eigvals,
                               const double lambda,
                               const double alpha) {
  /* ---------------------------------------------------------------------------
   Helper function that shrinks the eigenvalues.
   Shrinkage is that of the rotation equivariant alternative ridge estimator.
   - eigvals > A vector of containing the eigenvalues
   - lambda  > A double giving the ridge penalty
   - alpha   > A double giving the diagonal value in the rotaional equivariant
               target.
  --------------------------------------------------------------------------- */

  arma::vec seigvals = eigvals - lambda*alpha;
  return sqrt(lambda + 0.25f*pow(seigvals, 2.0f)) + 0.5f*seigvals;
}



arma::mat armaRidgePAnyTarget_OLD(const arma::mat & S,
                                  const arma::mat & target,
                                  const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general targets.
   Old ridge estimator using matrix inversion.
   - S      > A sample covariance matrix
   - target > The target matrix with the same size as S
   - lambda > The the ridge penalty
  --------------------------------------------------------------------------- */

  const int n = S.n_cols;
  const arma::mat E = symmatl(S) - lambda * symmatl(target);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, 0.25f*E*E + lambda*arma::eye(n, n), "dc");

  return inv_sympd(0.5f*E + eigvec*diagmat(sqrt(eigval))*eigvec.t());
}



// [[Rcpp::export(.armaRidgePAnyTarget)]]
arma::mat armaRidgePAnyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general targets.
   Not using matrix inversion.
   - S      > A sample covariance matrix
   - target > The target matrix with the same size as S
   - lambda > The ridge penalty
   NOTE -- S should not contain any NAs, Infs, or NaNs!
  --------------------------------------------------------------------------- */

  if (lambda == arma::datum::inf) {  //
    return target;
  }

  const int n = S.n_cols;
  const double inv_lambda = 1.0f/lambda;
  const arma::mat E = S - lambda * target;

  arma::vec eigval;
  arma::mat eigvec = E*E;  // Temporarily set to E*E

  // Return target if E*E contain infinite values due to large lambda
  // Usually happens for lambda greater than or on the order of 1e154
  if (!eigvec.is_finite()) {
    return target;
  }

  // The following line overwrite eigvec
  arma::eig_sym(eigval, eigvec, 0.25f*eigvec + lambda*arma::eye(n, n), "dc");
  if (any(eigval < 0)) { // Throw error if matrix square root is not defined.
    Rcpp::stop("Eigenvalues are not all positive. lambda is too small.");
  }

  eigval = inv_lambda*sqrt(eigval);  // Take the matrix square root and scale

  return symmatl((-0.5f*inv_lambda)*E + eigvec*diagmat(eigval)*eigvec.t());
}



arma::mat armaRidgePRotationInvariantTarget_OLD(const arma::mat & S,
                                                const double alpha,
                                                const double lambda) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for rotational equivariant targets.
   Old ridge estimator using matrix inversion.
   - S      > the sample covariance matrix
   - alpha  > The scaling of the identity matrix
   - lambda > The ridge penalty
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, symmatl(S), "dc");  // Eigen decomposition

  eigvals = armaEigShrink(eigvals, lambda, alpha);

  return inv_sympd(eigvecs*diagmat(eigvals)*eigvecs.t());
}



// [[Rcpp::export(.armaRidgePRotationInvariantTarget)]]
arma::mat armaRidgePRotationInvariantTarget(const arma::mat & S,
                                            const double alpha,
                                            const double lambda) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in the rotational equivariant case.
   Not using matrix inversion. When the shrunken eigen-values become infinite
   due to floating point errors, the function returns the target matrix.
   Usually happens for lambda >= 1e154
   - S      > A sample covariance matrix
   - alpha  > The scaling of the identity matrix.
   - lambda > The ridge penalty. Can be set to Inf (on the R side)
  --------------------------------------------------------------------------- */

  if (lambda == arma::datum::inf) {
    const int p = S.n_rows;
    return alpha*arma::eye<arma::mat>(p, p);
  }

  if (lambda <= 0) {
    Rcpp::stop("The penalty (lambda) must be strictly postive");
  }

  const double inv_lambda = 1.0f/lambda;
  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, S, "dc");

  // Eigenvalue shrinkage + addition to avoid inversion
  eigvals = inv_lambda*(armaEigShrink(eigvals,lambda,alpha) - eigvals) + alpha;

  // Throw error any eigenvalues are negative.
  if (any(eigvals < 0)) {
    Rcpp::stop("Eigenvalues are not all positive. lambda is too small.");
  }

  // Return target if shrunken evals are infinite and lambda is "large"
  // Usually happens for lambda >= 1e154
  if (!eigvals.is_finite() && lambda > 1) {
    const int p = S.n_rows;
    return alpha*arma::eye<arma::mat>(p, p);
  }

  // Transform back and return results
  return symmatl(eigvecs*diagmat(eigvals)*eigvecs.t());
}



// [[Rcpp::export(.armaRidgeP)]]
arma::mat armaRidgeP(const arma::mat & S,
                     const arma::mat & target,
                     const double lambda) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in C++.
   - S      > The sample covariance matrix (a numeric matrix on the R side)
   - target > Target matrix (a numeric matrix on the R side, same size as S)
   - lambda > The penalty (a numeric of length one on the R side)
  --------------------------------------------------------------------------- */

  const int n = S.n_rows;
  const double alpha = target(0, 0);
  const arma::mat alphaI = alpha*arma::eye<arma::mat>(n, n);
  if (arma::all(arma::all(target == alphaI))) {
    return armaRidgePRotationInvariantTarget(S, alpha, lambda);
  } else {
    return armaRidgePAnyTarget(S, target, lambda);
  }
}



////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

  FUSED RIDGE ESTIMATOR

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export(.armaFusedUpdateI)]]
arma::mat armaFusedUpdateI(int g0,
                           const Rcpp::List & Plist,
                           const Rcpp::List & Slist,
                           const Rcpp::List & Tlist,
                           const arma::vec & ns,
                           const double & lambda,
                           arma::mat & lambdaF) {
  /* ---------------------------------------------------------------------------
   "Update" the g0'th covariance matrix according to the first fusion
   update scheme and get the fused ridge estimate hereof.
   - g0          > An integer giving the class estimate to be updated.
   - Plist       > A list of length G of matrices giving the current precision
                   estimates.
   - Slist       > A list of length G of sample correlation matrices the same
                   size as those of Plist.
   - Tlist       > A list of length G of target matrices the same size
                   as those of Plist
   - ns          > A vector of length G giving the sample sizes.
   - lambda      > The ridge penalty (a postive number).
   - lambdaF     > A G by G symmetric adjacency matrix giving the fused
                   penalty graph with non-negative entries where
                   lambdaF[g1, g2] determine the (rate of) shrinkage
                   between estimates in classes corresponding to Slist[g1]
                   and Slist[g1].

   NOTE1: The C++ implmentation of .fusedUpdateI.
   NOTE2: The index starts from zero in C++. This is exported to the R side!
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double a = (sum(lambdaF.row(g0)) + lambda)/(ns[g0]);

  arma::mat Sbar = Rcpp::as<arma::mat>(Rcpp::wrap(Slist[g0]));
  arma::mat Tbar = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g0]));
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     arma::mat O = Rcpp::as<arma::mat>(Rcpp::wrap(Plist[g]));
     arma::mat T = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g]));
     Sbar -= (lambdaF(g, g0)/ns(g0))*(O - T);
  }
  return armaRidgeP(Sbar, Tbar, a);
}



arma::mat armaFusedUpdateIC(int g0,
                            const arma::cube & Pcube,
                            const arma::cube & Scube,
                            const arma::cube & Tcube,
                            const arma::vec ns,
                            const double lambda,
                            arma::mat lambdaF) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateI above with arma::cube instead of Rcpp::List in the
   arguments.
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double a = (sum(lambdaF.row(g0)) + lambda)/ns[g0];

  arma::mat Sbar = Scube.slice(g0);
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Sbar -= (lambdaF(g, g0)/ns(g0))*(Pcube.slice(g) - Tcube.slice(g));
  }
  return armaRidgeP(Sbar, Tcube.slice(g0), a);
}



// [[Rcpp::export(.armaFusedUpdateII)]]
arma::mat armaFusedUpdateII(int g0,
                            const Rcpp::List & Plist,
                            const Rcpp::List & Slist,
                            const Rcpp::List & Tlist,
                            const arma::vec ns,
                            const double lambda,
                            arma::mat lambdaF) {
  /* ---------------------------------------------------------------------------
   "Update" the g0'th covariance matrix according to the second fusion update
   scheme use the regular ridge estimate hereof.
   - g0          > An integer giving the class estimate to be updated.
   - Plist       > A list of length G of matrices giving the current precision
                   estimates.
   - Slist       > A list of length G of sample correlation matrices the same
                   size as those of Plist.
   - Tlist       > A list of length G of target matrices the same size
                   as those of Plist
   - ns          > A vector of length G giving the sample sizes.
   - lambda      > The ridge penalty (a postive number).
   - lambdaF  > A G by G symmetric adjacency matrix giving the fused
                   penalty graph with non-negative entries where
                   lambdaF[g1, g2] determine the (rate of) shrinkage
                   between estimates in classes corresponding to Slist[g1]
                   and Slist[g1].
    NOTE1: The C++ implmentation of .fusedUpdateII.
    NOTE2: The indexing starts in zero.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double rowsum = sum(lambdaF.row(g0));
  const double a = (rowsum + lambda)/ns[g0];
  const double b = (rowsum + lambda - 1)/ns(g0);

  arma::mat Sbar = Rcpp::as<arma::mat>(Rcpp::wrap(Slist[g0]));
  arma::mat Tbar = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g0]));
  const int p = Sbar.n_rows;

  arma::mat Psum = arma::zeros(p, p);
  arma::mat Tsum = Psum;
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     arma::mat P = Rcpp::as<arma::mat>(Rcpp::wrap(Plist[g]));
     arma::mat T = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g]));
     Psum += lambdaF(g0, g)*P;
     Tsum += (lambdaF(g0, g)/ns(g0))*T;
  }

  Sbar += b*Psum + Tsum;
  Tbar += Psum;
  return armaRidgeP(Sbar, Tbar, a);
}



arma::mat armaFusedUpdateIIC(int g0,
                             const arma::cube & Pcube,
                             const arma::cube & Scube,
                             const arma::cube & Tcube,
                             const arma::vec ns,
                             const double lambda,
                             arma::mat lambdaF) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateII above with arma::cube instead of Rcpp::List in the
   arguments.
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  const int p = Scube.n_rows;

  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double rowsum = sum(lambdaF.row(g0));
  const double a = (rowsum + lambda)/ns[g0];
  const double b = (rowsum + lambda - 1)/ns(g0);

  arma::mat Sbar = Scube.slice(g0);
  arma::mat Tbar = Tcube.slice(g0);
  arma::mat Psum = arma::zeros(p, p);
  arma::mat Tsum = Psum;
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Psum += lambdaF(g0, g)*Pcube.slice(g);
     Tsum += (lambdaF(g0, g)/ns(g0))*Tcube.slice(g);
  }

  Sbar += b*Psum + Tsum;
  Tbar += Psum;
  return armaRidgeP(Sbar, Tbar, a);
}



// [[Rcpp::export(.armaFusedUpdateIII)]]
arma::mat armaFusedUpdateIII(int g0,
                             const Rcpp::List & Plist,
                             const Rcpp::List & Slist,
                             const Rcpp::List & Tlist,
                             const arma::vec & ns,
                             const double & lambda,
                             arma::mat & lambdaF) {
  /* ---------------------------------------------------------------------------
   "Update" the g0'th covariance matrix according to the third fusion update
   scheme use the regular ridge estimate hereof.
   - g0          > An integer giving the class estimate to be updated.
   - Plist       > A list of length G of matrices giving the current precision
                   estimates.
   - Slist       > A list of length G of sample correlation matrices the same
                   size as those of Plist.
   - Tlist       > A list of length G of target matrices the same size
                   as those of Plist
   - ns          > A vector of length G giving the sample sizes.
   - lambda      > The ridge penalty (a postive number).
   - lambdaF  > A G by G symmetric adjacency matrix giving the fused
                   penalty graph with non-negative entries where
                   lambdaF[g1, g2] determine the (rate of) shrinkage
                   between estimates in classes corresponding to Slist[g1]
                   and Slist[g1].
    NOTE: The C++ implementaiton of .fusedUpdateIII.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double lambdasum = sum(lambdaF.row(g0)) + lambda;
  const double a = lambdasum/ns[g0];

  arma::mat Sbar = Rcpp::as<arma::mat>(Rcpp::wrap(Slist[g0]));
  arma::mat Tbar = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g0]));
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     arma::mat P = Rcpp::as<arma::mat>(Rcpp::wrap(Plist[g]));
     arma::mat T = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g]));
     Tbar += (lambdaF(g0, g)/lambdasum)*(P - T);
  }

  return armaRidgeP(Sbar, Tbar, a);
}



arma::mat armaFusedUpdateIIIC(int g0,
                              const arma::cube & Pcube,
                              const arma::cube & Scube,
                              const arma::cube & Tcube,
                              const arma::vec ns,
                              const double lambda,
                              arma::mat lambdaF) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateIII with arma::cube instead of Rcpp::List
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  lambdaF(g0, g0) = 0;  // Make sure entry (g0, g0) is zero
  const double lambdasum = sum(lambdaF.row(g0)) + lambda;
  const double a = lambdasum/ns[g0];
  arma::mat Tbar = Tcube.slice(g0);
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Tbar += (lambdaF(g0, g)/lambdasum)*(Pcube.slice(g) - Tcube.slice(g));
  }

  return armaRidgeP(Scube.slice(g0), Tbar, a);
}



// [[Rcpp::export(.armaRidgeP_fused)]]
Rcpp::List armaRidgeP_fused(const Rcpp::List & Slist,
                            const arma::vec & ns,
                            const Rcpp::List & Tlist,
                            const double & lambda,
                            arma::mat & lambdaF,
                            const Rcpp::List & Plist,
                            const int maxit = 100,
                            const double eps = 1e-5,
                            const bool verbose = false) {
  /* ---------------------------------------------------------------------------
   The fused ridge estimate workhorse function for a given lambda and
   lambdaF.
   - Slist   > A list of length G of sample correlation matrices the same size
               as those of Plist.
   - Tlist   > A list of length G of target matrices the same size
               as those of Plist. Default is given by default.target.
   - ns      > A vector of length G giving the sample sizes.
   - lambda  > The ridge penalty (a stictly postive number)
   - lambdaF > The G by G symmetric adjacency matrix fused penalty graph
               with non-negative entries where lambdaF[g1, g2] determine the
               retainment of similarities between estimates in classes
               corresponding to Slist[g1] and Slist[g1].
   - Plist   > A list of length G of symmetric p.d. matrices serves as initial
               estimates of the algorithm.
   - maxit   > integer. The maximum number of interations, default is 100.
   - eps     > numeric. A positive convergence criterion. Default is 1e-5.
   - verbose > logical. Should the function print extra info. Defaults to false.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();

  // Initialize
  const double lambdasize = lambda + accu(lambdaF);
  double delta;
  arma::vec diffs = arma::ones(G);  // Vector of ones, will be overwritten
  arma::mat tmp;
  Rcpp::List Plist_out = Rcpp::clone(Plist);

  for (int i = 0; i < maxit; ++i) {
    for (int g = 0; g < G; ++g) {
      tmp = Rcpp::as<arma::mat>(Plist_out(g));
      if (lambdasize < 1e50) {
        // Update I is faster but unstable for very large lambda
        Plist_out(g) = armaFusedUpdateI(g, Plist_out, Slist, Tlist, ns,
                                        lambda, lambdaF);
      } else {
        // Update III is slower but more stable for very large lambda
        Plist_out(g) = armaFusedUpdateIII(g, Plist_out, Slist, Tlist, ns,
                                          lambda, lambdaF);
      }
      diffs(g) = pow(norm(Rcpp::as<arma::mat>(Plist_out(g)) - tmp, "fro"), 2.0);
    }
    delta = max(diffs);
    if (delta > eps) {
      if (verbose) {
        Rprintf("i = %-3d | max diffs = %0.10f\n", i + 1, delta);
      }
    } else {
      if (verbose) {
        Rprintf("Converged in %d iterations.\n", i + 1);
      }
      break;
    }
  }
  return Plist_out;
}



////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

   GENERAL SIMULATION TOOLS

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
arma::mat rmvnormal(const int n, arma::rowvec mu, arma::mat sigma) {
  /* ---------------------------------------------------------------------------
   Simulate from multivariate normal distribution with mean mu and covariance
   sigma. Returns a matrix of size n by length(mu) of observations.
   - n     > An integer giving the number of samples to simulate.
   - mu    > A vector giving the population mean.
   - sigma > A matrix giving the population covariance matrix.
   (Copied taken from package GMCM)
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



arma::mat armaRWishartSingle(const arma::mat L,
                             const double nu,
                             const int p) {
  /* ---------------------------------------------------------------------------
    Wishart simulation. Simulate a single Wishart distributed matrix
    (Taken from package AEBilgrau/correlateR)
    L  > A cholesky decomposition matrix
    nu > The degrees of freedom
    p  > The dimension of the distribution
    (Copied from package AEBilgrau/correlateR)
  --------------------------------------------------------------------------- */

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



// [[Rcpp::export(.armaRWishart)]]
arma::cube armaRWishart(const int n,
                        const arma::mat & sigma,
                        const double nu) {
  /* ---------------------------------------------------------------------------
    Simulate n wishart distributed matrices.
    n     > A integer giving the number of matrices to generate
    sigma > The scale matrix in the Wishart distribution
    nu    > The degrees of freedom in the Wishart distribution
    (Copied from package AEBilgrau/correlateR)
  --------------------------------------------------------------------------- */

  const int p = sigma.n_cols;
  const arma::mat L = arma::chol(sigma);
  arma::cube ans(p, p, n);
  for (int g = 0; g < n; ++g) {
    ans.slice(g) = armaRWishartSingle(L, nu, p);
  }
  return ans;
}



// [[Rcpp::export(.armaRInvWishart)]]
arma::cube armaRInvWishart(const int n,
                           const arma::mat & psi,
                           const double nu) {
  /* ---------------------------------------------------------------------------
    Simulate n inverse Wishart distributed matrices
    n     > A integer giving the number of matrices to generate
    psi   > The scale matrix in the inverse Wishart distribution
    nu    > The degrees of freedom in the inverse Wishart distribution
    (Copied from package AEBilgrau/correlateR)
  --------------------------------------------------------------------------- */

  const int p = psi.n_cols ;
  const arma::mat L = arma::chol(inv(psi));
  arma::cube ans(p, p, n);
  for (int g = 0; g < n; ++g) {
    ans.slice(g) = inv(armaRWishartSingle(L, nu, p));
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
               B1 <- armaRidgePAnyTarget(S, target, lambda),
               C1 <- armaRidgeP(S, target, lambda),
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
               B2 <- armaRidgePRotationInvariantTarget(S, 0.0, lambda),
               C2 <- armaRidgeP(S, target, lambda))
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
               B3 <- armaRidgePRotationInvariantTarget(S, target[1,1], lambda),
               C3 <- armaRidgeP(S, target, lambda))
stopifnot(all.equal(unname(A3), unname(B3)))
stopifnot(all.equal(unname(A3), unname(C3)))



library("rags2ridges")
library("microbenchmark")
S  <- createS(n = 10, p = 500)
target <- default.target(S, type = "DEPV")#default.target(S)#
lambda <- 2
system.time(B <- ridgeSArma(S, target = target, lambda = lambda))
microbenchmark(A <- ridgeS(S, target = target, lambda = lambda),
               B <- armaRidgeP(S, target = target, lambda = lambda),
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
Plist <- createS(n = c(100, 100), p = 10)
Slist <- createS(n = c(100, 100), p = 10)
Tlist <- default.target.fused(Slist, ns)
g0 <- 1

res <- microbenchmark(
  A = armaFusedUpdate2(g0, Plist, Slist, Tlist, ns, 1, matrix(1, 2, 2)),
  B = armaFusedUpdate(g0, Plist, Slist, Tlist, ns, 1, matrix(1, 2, 2)),
  C = rags2ridges:::.fusedUpdate(g0, Plist, Slist, Tlist, ns, 1, matrix(1, 2, 2)),
  times = 10000)
boxplot(res)

all.equal(A, C)
all.equal(A2, A)
all.equal(A, B)
*/

