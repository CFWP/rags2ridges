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

// [[Rcpp::export(NLL)]]
double NLL(const arma::mat S, const arma::mat P) {
  /* The negative loglikelihood */
  double logdet;
  double sign;
  log_det(logdet, sign, P);
  if (sign < 0) {
    Rcpp::warning("Supplied precision matrix is not postive definite.");
    return arma::datum::inf;
  }
  return -logdet + accu(S % P);
}



// [[Rcpp::export(PNLL)]]
double PNLL(const arma::mat S, const arma::mat P, const arma::mat T,
            const double lambda) {
  /* The penalized negative loglikelihood */
  return NLL(S, P) + 0.5*lambda*pow(arma::norm(P - T, "fro"), 2.0);
}



// [[Rcpp::export(NLL.fused)]]
double NLL_fused(const Rcpp::List Slist, const Rcpp::List Plist,
                const arma::vec ns) {
  /* ---------------------------------------------------------------------------
   Function that computes the value of the (negative) combined log-likelihood
   - Slist > A list sample covariance matrices for each class
   - Plist > A list of the same length as (Slist) of precision matrices
             (possibly regularized inverse covariance or correlation matrices)
   - ns    > A vector of sample sizes of the same length as Slist.
  --------------------------------------------------------------------------- */

  const int G = ns.size();
  double nll = 0;
  for (int i = 0; i < G; ++i) {
    arma::mat Si = Slist[i];
    arma::mat Pi = Plist[i];
    nll += ns[i]*NLL(Si, Pi);
  }
  return nll;
}



// [[Rcpp::export(PNLL.fused)]]
double PNLL_fused(const Rcpp::List Slist, const Rcpp::List Plist,
                  const arma::vec ns, const Rcpp::List Tlist,
                  const arma::mat lambda) {
  /* ---------------------------------------------------------------------------
   Function that computes the value of the penalized (negative) fused
   log-likelihood.
   - Slist  > A list sample covariance matrices for each class
   - Plist  > A list of the same length as (Slist) of precision matrices
              (possibly regularized inverse covariance or correlation matrices)
   - Tlist  > A list of target matrices
   - ns     > A vector of sample sizes of the same length as Slist.
   - lambda > The penalty matrix
  --------------------------------------------------------------------------- */

  const int G = ns.size();
  double pnll = NLL_fused(Slist, Plist, ns);
  for (int i = 0; i < G; i++) {
    arma::mat Pi = Plist[i];
    arma::mat Ti = Tlist[i];
    arma::mat diff = Pi - Ti;
    pnll += 0.5*lambda(i,i)*pow(arma::norm(diff, "fro"), 2.0);
    for (int j = 0; j < i; j++) {
      arma::mat Pj = Plist[j];
      arma::mat Tj = Tlist[j];
      pnll += 0.25*lambda(i,j)*pow(arma::norm(diff - Pj + Tj, "fro"), 2.0);
    }
  }
  return pnll;
}



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



// [[Rcpp::export(.armaPooledP)]]
arma::mat armaPooledP(const Rcpp::List & Plist,  // List of precision matrices
                      const Rcpp::NumericVector ns,
                      const int mle = 0) {
  /* ---------------------------------------------------------------------------
   Function to compute the pooled precision estimate. Returns a numeric matrix
   giving the pooled precision matrix.
   Simply a wrapper for armaPooledS and matrix inversion.
   - Plist > A list of (perhaps estimated) precision matrices.
   - nu    > A numeric vector giving the number of samples corresponding
             to each scatter matrix.
   - mle   > A logical (or integer) of length one equalling FALSE (0) or
             TRUE (1). If FALSE, the bias corrected ML estimate is used.
             If TRUE the ML estimate is used.
  --------------------------------------------------------------------------- */

  const int G = Plist.size();
  const int imle = 1 - mle;
  const double rdenum = 1.0f/(sum(ns) - G * imle);
  arma::mat S0 = Plist[0];
  S0 = (ns[0] - imle)*arma::inv_sympd(S0);
  for (int i = 1; i < G; ++i) {  // Loop through remaning elements. Note i = 1.
    arma::mat Pi = Plist[i];
    S0 += (ns[i] - imle)*arma::inv_sympd(Pi);
  }
  return arma::inv_sympd(rdenum*S0);
}



// [[Rcpp::export(.armaEigShrink)]]
arma::vec armaEigShrink(const arma::vec dVec,
                        const double lambda,
                        const double cons = 0) {
  /* ---------------------------------------------------------------------------
   - Function that shrinks the eigenvalues
   - Shrinkage is that of the rotation equivariant alternative ridge estimator
   - Main use is in avoiding expensive matrix square root when choosing a
     target that leads to a rotation equivariant version of the alternative
     ridge estimator
   - dVec   > numeric vector containing the eigenvalues of a matrix S
   - lambda > penalty parameter
   - const  > a constant, default = 0
   --------------------------------------------------------------------------- */

  arma::vec Evector = 0.5 * (dVec - lambda * cons);
  return sqrt(lambda + pow(Evector, 2.0)) + Evector;
}



// [[Rcpp::export(.armaEigShrinkAnyTarget)]]
arma::vec armaEigShrinkAnyTarget(const arma::mat & S,
                                 const arma::mat & target,
                                 const double lambda) {
  /* ---------------------------------------------------------------------------
  - Function that shrinks the eigenvalues
  - Shrinkage is that of the alternative ridge estimator under a general target
  - Main use is in avoiding expensive Schur-approach to computing the
    matrix square root
  - S      > A sample covariance matrix
  - target > Target matrix of same dimensions as S
  - lambda > penalty parameter
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs = S - lambda * target;
  eig_sym(eigvals, eigvecs, eigvecs, "dc");
  eigvals = 0.5 * eigvals;
  arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));
  return (sqroot + eigvals);
}



// [[Rcpp::export(.armaEigShrinkArchI)]]
arma::vec armaEigShrinkArchI(const arma::vec dVec,
                             const double lambda,
                             const double cons) {
  /* ---------------------------------------------------------------------------
   - Function that shrinks the eigenvalues
   - Shrinkage is that of the rotation equivariant Archetypal I estimator
   - dVec   > numeric vector containing the eigenvalues of a matrix S
   - lambda > penalty parameter
   - const  > a constant
   --------------------------------------------------------------------------- */

  arma::vec Evector = (1 - lambda) * dVec + (lambda * (1.0/cons));
  return Evector;
}




////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

  TOOLS FOR THE REGULAR (NON-FUSED) RIDGE ESTIMATOR OF THE CORE MODULE

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////

inline arma::mat rev_eig(const arma::vec eigval, const arma::mat eigvec) {
  /* ---------------------------------------------------------------------------
   "Reverse" the eigen decomposition, i.e. perform the multiplcation
   - eigval > a vector of eigenvalues
   - eigvec > a matrix of corresponding eigenvectors
  --------------------------------------------------------------------------- */

  return eigvec*diagmat(eigval)*eigvec.t();  // Could be more efficient
}



// [[Rcpp::export(.armaRidgePAnyTarget)]]
arma::mat armaRidgePAnyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda,
                              int invert = 2) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general/arbitrary targets.
   Depending on the value of "invert"" using matrix inversion (via
   diagonalization) or avoiding it.
   - S      > A sample covariance matrix. Should not contain NAs, Infs, or NaNs!
   - target > The target matrix with the same size as S
   - lambda > The the ridge penalty
   - invert > integer. Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic" (default).
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs = S - lambda*target;
  if (!eigvecs.is_finite()) {
    return target;
  }
  eig_sym(eigvals, eigvecs, eigvecs, "dc");
  eigvals = 0.5*eigvals;
  arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));

  // Return target if shrunken evals are infinite and lambda is "large"
  // Usually happens for lambda >= 1e154
  if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
    return target;
  }

  arma::vec D_inv = 1.0/(sqroot + eigvals); // inversion diagonal

  if (invert == 2) {   // Determine to invert or not
    if (lambda > 1) {  // Generally, don't use inversion for "large" lambda
      invert = 0;
    } else {
      if (!D_inv.is_finite()) {
        invert = 0;
      } else {
        invert = 1;
      }
    }
  }

  // Determine to invert or not
  if (invert == 1) {
    return rev_eig(D_inv, eigvecs);  // Proper inversion
  } else {
    arma::vec D_noinv = (sqroot - eigvals)/lambda; // inversion-less diagonal
    return rev_eig(D_noinv, eigvecs);  // Inversion by proposion
  }

}



// [[Rcpp::export(.armaRidgePScalarTarget)]]
arma::mat armaRidgePScalarTarget(const arma::mat & S,
                                 const double alpha,
                                 const double lambda,
                                 int invert = 2) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for rotational equivariant targets.
   Depending on the value of "invert"" using matrix inversion (via
   diagonalization) or avoiding it.
   - S      > A sample covariance matrix.
   - alpha  > The scaling of the identity matrix. Shoud not contain NaNs, Infs,
              or NA.s
   - lambda > The ridge penalty. Can be set to Inf (on the R side)
   - invert > Should the estimate be computed using inversion?
              0 = "no", 1 = "yes", 2 = "automatic", (default).
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, S, "dc");

  eigvals = 0.5*(eigvals - lambda*alpha);
  arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));

  // Return target if shrunken evals are infinite and lambda is "large"
  // Usually happens for lambda >= 1e154
  if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
    const int p = S.n_rows;
    return alpha*arma::eye<arma::mat>(p, p);
  }

  arma::vec D_inv = 1.0/(sqroot + eigvals); // inversion diagonal

  if (invert == 2) {   // Determine to invert or not
    if (lambda > 1) {  // Generally, don't use inversion for "large" lambda
      invert = 0;
    } else {
      if (!D_inv.is_finite()) {
        invert = 0;
      } else {
        invert = 1;
      }
    }
  }

  // Determine to invert or not
  if (invert == 1) {
    return rev_eig(D_inv, eigvecs);  // Proper inversion
  } else {
    arma::vec D_noinv = (sqroot - eigvals)/lambda; // inversion-less diagonal
    return rev_eig(D_noinv, eigvecs);  // Inversion by proposition
  }

}



// [[Rcpp::export(.armaRidgeP)]]
arma::mat armaRidgeP(const arma::mat & S,
                     const arma::mat & target,
                     const double lambda,
                     int invert = 2) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in C++. Wrapper for the subroutines
   - S      > The sample covariance matrix (a numeric matrix on the R side)
   - target > Target matrix (a numeric matrix on the R side, same size as S)
   - lambda > The penalty (a numeric of length one on the R side)
   - invert > Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic", (default).
  --------------------------------------------------------------------------- */

 if (lambda <= 0) {
    Rcpp::stop("The penalty (lambda) must be strictly postive");
  }

  if (lambda == arma::datum::inf) {
    return target;
  }

  const int p = S.n_rows;
  const double alpha = target(0, 0);
  const arma::mat alphaI = alpha*arma::eye<arma::mat>(p, p);

  if (arma::all(arma::all(target == alphaI))) {
    return armaRidgePScalarTarget(S, alpha, lambda, invert);
  } else {
    return armaRidgePAnyTarget(S, target, lambda, invert);
  }

}




////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

 TOOLS FOR THE FUSED RIDGE ESTIMATOR OF THE FUSED MODULE

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(.armaFusedUpdateI)]]
arma::mat armaFusedUpdateI(int g0,
                           const Rcpp::List & Plist,
                           const Rcpp::List & Slist,
                           const Rcpp::List & Tlist,
                           const arma::vec & ns,
                           const arma::mat & lambda) {
  /* ---------------------------------------------------------------------------
   Updated the g0'th precision estimate according to the first fused ridge
   update scheme.
   NOTE: The index starts from zero in C++. This is exported to the R side!
   NOTE: No check are made that the arguments are of correct type and format.
   - g0     > An integer between 0 and G giving the class estimate to be
              updated.
   - Plist  > A list of G matrices giving the current precision estimates.
   - Slist  > A list of G sample covariance matrices the same size as in Plist.
   - Tlist  > A list of length G of target matrices the same size as in Plist
   - ns     > A vector of length G giving the sample sizes.
   - lambda > A G by G non-negative symmetric adjacency penalty matrix. The
              diagonal entries are the ridge penalties and should be strictly
              postive. The off-diagonal entries are the non-negative fusion
              penalties. E.g. lambda[g1, g2] determine the (rate of) shrinkage
              between estimates in classes corresponding to g1 and g2.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  const double lambda_a = sum(lambda.row(g0))/ns[g0];

  arma::mat Sbar = Rcpp::as<arma::mat>(Rcpp::wrap(Slist[g0]));
  arma::mat Tbar = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g0]));
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     arma::mat O = Rcpp::as<arma::mat>(Rcpp::wrap(Plist[g]));
     arma::mat T = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g]));
     Sbar -= (lambda(g, g0)/ns[g0])*(O - T);
  }
  return armaRidgeP(Sbar, Tbar, lambda_a);
}



arma::mat armaFusedUpdateIC(int g0,
                            const arma::cube & Pcube,
                            const arma::cube & Scube,
                            const arma::cube & Tcube,
                            const arma::vec ns,
                            const arma::mat lambda) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateI above with arma::cube instead of Rcpp::List in the
   arguments.
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  const double lambda_a = sum(lambda.row(g0))/ns[g0];

  arma::mat Sbar = Scube.slice(g0);
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Sbar -= (lambda(g, g0)/ns[g0])*(Pcube.slice(g) - Tcube.slice(g));
  }
  return armaRidgeP(Sbar, Tcube.slice(g0), lambda_a);
}



// [[Rcpp::export(.armaFusedUpdateII)]]
arma::mat armaFusedUpdateII(int g0,
                            const Rcpp::List & Plist,
                            const Rcpp::List & Slist,
                            const Rcpp::List & Tlist,
                            const arma::vec ns,
                            const arma::mat lambda) {
  /* ---------------------------------------------------------------------------
   Updates the g0'th precision estimate according to the second fusion update
   scheme.
   NOTE: The index starts from zero in C++. This is exported to the R side!
   - g0     > An integer between 0 and G giving the class to be updated.
   - Plist  > A list of G matrices giving the current precision estimates.
   - Slist  > A list of G sample covariance matrices the same size as in Plist.
   - Tlist  > A list of length G of target matrices the same size as in Plist
   - ns     > A vector of length G giving the sample sizes.
   - lambda > A G by G non-negative symmetric adjacency penalty matrix. The
              diagonal entries are the ridge penalties and should be strictly
              postive. The off-diagonal entries are the non-negative fusion
              penalties. E.g. lambda[g1, g2] determine the (rate of) shrinkage
              between estimates in classes corresponding to g1 and g2.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  const double lambdasum = sum(lambda.row(g0));
  const double a = lambdasum/ns[g0];
  const double b = (lambdasum - 1.0)/ns[g0];

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
     Psum += lambda(g0, g)*P;
     Tsum += (lambda(g0, g)/ns[g0])*T;
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
                             const arma::mat lambda) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateII above with arma::cube instead of Rcpp::List in the
   arguments.
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  const int p = Scube.n_rows;

  const double lambdasum = sum(lambda.row(g0));
  const double a = lambdasum/ns[g0];
  const double b = (lambdasum - 1.0)/ns(g0);

  arma::mat Sbar = Scube.slice(g0);
  arma::mat Tbar = Tcube.slice(g0);
  arma::mat Psum = arma::zeros(p, p);
  arma::mat Tsum = Psum;
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Psum += lambda(g0, g)*Pcube.slice(g);
     Tsum += (lambda(g0, g)/ns[g0])*Tcube.slice(g);
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
                             const arma::mat & lambda) {
  /* ---------------------------------------------------------------------------
   Updates the g0'th precision matrix according to the third fusion update
   scheme.
   NOTE: The index starts from zero in C++. This is exported to the R side!
   - g0     > An integer giving the class estimate to be updated.
   - Plist  > A list of G matrices giving the current precision estimates.
   - Slist  > A list of G sample covariance matrices the same size as in Plist.
   - Tlist  > A list of length G of target matrices the same size as in Plist
   - ns     > A vector of length G giving the sample sizes.
   - lambda > A G by G non-negative symmetric adjacency penalty matrix. The
              diagonal entries are the ridge penalties and should be strictly
              postive. The off-diagonal entries are the non-negative fusion
              penalties. E.g. lambda[g1, g2] determine the (rate of) shrinkage
              between estimates in classes corresponding to g1 and g2.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  const double a = sum(lambda.row(g0));

  arma::mat Sbar = Rcpp::as<arma::mat>(Rcpp::wrap(Slist[g0]));
  arma::mat Tbar = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g0]));
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }

     double factor = arma::is_finite(lambda(g0, g)) ? lambda(g0, g)/a : 1;

     arma::mat P = Rcpp::as<arma::mat>(Rcpp::wrap(Plist[g]));
     arma::mat T = Rcpp::as<arma::mat>(Rcpp::wrap(Tlist[g]));
     Tbar += factor*(P - T);
  }

  return armaRidgeP(Sbar, Tbar, a/ns[g0]);
}



arma::mat armaFusedUpdateIIIC(int g0,
                              const arma::cube & Pcube,
                              const arma::cube & Scube,
                              const arma::cube & Tcube,
                              const arma::vec & ns,
                              const arma::mat & lambda) {
  /* ---------------------------------------------------------------------------
   As armaFusedUpdateIII with arma::cube instead of Rcpp::List
  --------------------------------------------------------------------------- */

  const int G = Scube.n_slices;
  const double lambdasum = sum(lambda.row(g0));
  const double a = lambdasum/ns[g0];
  arma::mat Tbar = Tcube.slice(g0);
  for (int g = 0; g < G; ++g) {
     if (g == g0) {
       continue;
     }
     Tbar += (lambda(g0, g)/lambdasum)*(Pcube.slice(g) - Tcube.slice(g));
  }

  return armaRidgeP(Scube.slice(g0), Tbar, a);
}



// [[Rcpp::export(.armaRidgeP.fused)]]
Rcpp::List armaRidgeP_fused(const Rcpp::List & Slist,
                            const arma::vec & ns,
                            const Rcpp::List & Tlist,
                            const arma::mat & lambda,
                            const Rcpp::List & Plist,
                            const int maxit = 100,
                            const double eps = 1e-7,
                            const bool relative = true,
                            const bool verbose = false) {
  /* ---------------------------------------------------------------------------
   The fused ridge estimate workhorse function for a given lambda.
   - Slist   > A list of G of sample covariance matrices of the same size.
   - Tlist   > A list of G of p.d. target matrices the same size as in Slist.
   - ns      > A vector of length G giving the sample sizes.
   - lambda  > A G by G non-negative symmetric adjacency penalty matrix. The
               diagonal entries are the ridge penalties and should be strictly
               postive. The off-diagonal entries are the non-negative fusion
               penalties. E.g. lambda[g1, g2] determine the (rate of) shrinkage
               between estimates in classes corresponding to g1 and g2.
   - Plist   > A list of length G of symmetric p.d. matrices that serves as
               initial estimates of the algorithm.
   - maxit   > integer. The maximum number of interations, default is 100.
   - eps     > numeric. A positive convergence criterion. Default is 1e-7.
   - relative > Devide
   - verbose > logical. Should the function print extra info. Defaults to false.
  --------------------------------------------------------------------------- */

  const int G = Slist.size();
  double delta;
  const arma::colvec lambda_colsum = sum(lambda, 1);  // row sums
  arma::vec diffs = arma::ones(G);  // Vector of ones, will be overwritten
  arma::mat tmp;
  Rcpp::List Plist_out = Rcpp::clone(Plist);

  for (int i = 0; i < maxit; ++i) {
    for (int g = 0; g < G; ++g) {

      tmp = Rcpp::as<arma::mat>(Plist_out(g));
      if (lambda_colsum[g] < 1) {
        // Update I is more stable for very large lambda
        Plist_out(g) = armaFusedUpdateI(g, Plist_out, Slist, Tlist, ns, lambda);
      } else {
        // Update III is more stable for very large lambda
        Plist_out(g) = armaFusedUpdateIII(g, Plist_out, Slist, Tlist,ns,lambda);
      }
      diffs(g) = pow(norm(Rcpp::as<arma::mat>(Plist_out(g)) - tmp, "fro"), 2.0);
      if (relative) {
        diffs(g) /=  pow(norm(Rcpp::as<arma::mat>(Plist_out(g)), "fro"), 2.0);
      }
    }
    delta = max(diffs);

    if (verbose) {
      Rprintf("i = %-3d | max diff = %-15.10e\n", i + 1, delta);
    }
    if (delta < eps) {
      if (verbose) {
        Rprintf("Converged in %d iterations, max diff < %1.2e.\n", i + 1, eps);
      }
      return Plist_out;
    }
  }

  Rcpp::warning("Max iterations (%d) hit.", maxit);
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
   NOTE: Copied from package GMCM
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
    L  > A cholesky decomposition matrix
    nu > The degrees of freedom
    p  > The dimension of the distribution
    NOTE: Copied from package AEBilgrau/correlateR
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
    NOTE: Copied from package AEBilgrau/correlateR
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
    NOTE: Copied from package AEBilgrau/correlateR
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
getwd()
source("./../../fnorm.R")
*/

