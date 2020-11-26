// We only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>

////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------

   AUXILIARY TOOLS

----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////

//' Evaluate the (penalized) (fused) likelihood
//'
//' Functions that evaluate the (penalized) (fused) likelihood.
//'
//' @aliases NLL PNLL NLL.fused PNLL.fused
//' @param S,Slist A (list of) positive semi definite sample covariance
//' matrices.
//' @param P,Plist A (list of) positive definite precision matrices.
//' @param T,Tlist A (list of) positive definite target matrices.
//' @param ns A \code{numeric} of sample sizes.
//' @param lambda A \code{numeric} penalty parameter. For the \code{.fused}
//' functions, this is a penalty \code{matrix}.
//' @return A single number.
//' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
//' Wessel N. van Wieringen
//' @seealso \code{\link{ridgeP}}, \code{\link{ridgeP.fused}}
//' @examples
//'
//' ns <- c(4,5)
//' Slist <- createS(n = ns, p = 5)
//' Plist <- list(diag(5), diag(2,5))
//' Tlist <- list(diag(5), diag(5))
//'
//' NLL(Slist[[1]], Plist[[1]])
//' PNLL(Slist[[1]], Plist[[1]], Tlist[[1]], lambda = 1)
//' NLL.fused(Slist, Plist, ns)
//' PNLL.fused(Slist, Plist, ns, Tlist, lambda = diag(2))
//'
//' @export
// [[Rcpp::export]]
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


//' @rdname NLL
//' @export
// [[Rcpp::export]]
double PNLL(const arma::mat S, const arma::mat P, const arma::mat T,
            const double lambda) {
  /* The penalized negative loglikelihood */
  return NLL(S, P) + 0.5*lambda*pow(arma::norm(P - T, "fro"), 2.0);
}

//' @rdname NLL
//' @export
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

//' @rdname NLL
//' @export
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


//' Core ridge precision estimators
//'
//' This is the interface to the \code{C++} implementations of the ridge
//' precision estimators. They perform core computations for many other
//' functions.
//'
//' @details
//' The functions are R-interfaces to low-level \code{C++} implementations
//' of the ridge estimators in the reference below
//' (cf. Lemma 1, Remark 6, Remark 7, and section 4.1 therein).
//'
//' \code{.armaRidgeP} is simply a wrapper (on the C++ side) for
//' \code{.armaRidgePAnyTarget} and \code{.armaRidgePScalarTarget} which are
//' the estimators for arbitrary and scalar targets, respectively.
//' The \code{invert} argument of the functions indicates if the computation
//' uses matrix inversion or not.
//'
//' Essentially, the functions compute
//'     \deqn{
//'       \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
//'         \left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
//'         \lambda_{a}\mathbf{T})^{2}\right]^{1/2} + \frac{1}{2}(\mathbf{S} -
//'         \lambda_{a}\mathbf{T})\right\}^{-1},
//'     }{%
//'       \Omega(\lambda) =
//'       \{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
//'         + 1/2(S - \lambda T)\}^{-1},
//'     }
//' if \code{invert == 1} or
//'     \deqn{
//'       \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
//'         \frac{1}{\lambda}\left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
//'                  \lambda_{a}\mathbf{T})^{2}\right]^{1/2} - \frac{1}{2}(\mathbf{S} -
//'         \lambda_{a}\mathbf{T})\right\}
//'     }{%
//'       \Omega(\lambda) =
//'       1/{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
//'         - 1/2(S - \lambda T)},
//'     }
//' if \code{invert == 0} using appropriate eigenvalue decompositions.
//' See the \R implementations in the example section below.
//'
//' @aliases .armaRidgeP .armaRidgePAnyTarget .armaRidgePScalarTarget armaRidgeP
//' armaRidgePAnyTarget armaRidgePScalarTarget
//' @param S A sample covariance \code{matrix}.
//' @param target A \code{numeric} symmetric positive definite target
//' \code{matrix} of the same size as \code{S}.
//' @param alpha The diagonal value on the scalar target matrix. A \code{double}
//' of length 1.
//' @param lambda The ridge penalty. A \code{double} of length 1.
//' @param invert An \code{integer}. Should the estimate be computed using
//' inversion?  Permitted values are \code{0L}, \code{1L}, and \code{2L} meaning
//' "no", "yes", and "automatic" (default).
//' @return Returns a symmetric positive definite \code{matrix} of the same size
//' as \code{S}.
//' @section Warning: The functions themselves perform no checks on the input.
//' Correct input should be ensured by wrappers.
//' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>
//' @seealso Used as a backbone in \code{\link{ridgeP}},
//' \code{\link{ridgeP.fused}}, etc.
//' @references van Wieringen, W.N. & Peeters, C.F.W. (2016).  Ridge Estimation
//' of Inverse Covariance Matrices from High-Dimensional Data, Computational
//' Statistics & Data Analysis, vol. 103: 284-303.  Also available as
//' arXiv:1403.0904v3 [stat.ME].
//' @keywords internal
//' @examples
//'
//' S <- createS(n = 3, p = 4)
//' tgt <- diag(4)
//' rags2ridges:::.armaRidgeP(S, tgt, 1.2)
//'
//' rags2ridges:::.armaRidgePAnyTarget(S, tgt, 1.2)
//' rags2ridges:::.armaRidgePScalarTarget(S, 1, 1.2)
//'
//'
//' ################################################################################
//' # The C++ estimators essentially amount to the following functions.
//' ################################################################################
//'
//' rev_eig <- function(evalues, evectors) { # "Reverse" eigen decomposition
//'   evectors %*% diag(evalues) %*% t(evectors)
//' }
//'
//' # R implementations
//'
//' # As armaRidgePScalarTarget Inv
//' rRidgePScalarTarget <- function(S, a, l, invert = 2) {
//'   ED <- eigen(S, symmetric = TRUE)
//'   eigvals <- 0.5*(ED$values - l*a)
//'   sqroot <- sqrt(l + eigvals^2)
//'
//'   if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
//'     return(diag(a, nrow(S)))
//'   }
//'
//'   D_inv <- 1.0/(sqroot + eigvals)
//'   D_noinv <- (sqroot - eigvals)/l
//'
//'   if (invert == 2) {   # Determine to invert or not
//'     if (l > 1) {  # Generally, use inversion for "small" lambda
//'       invert = 0;
//'     } else {
//'       invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
//'     }
//'   }
//'
//'   if (invert) {
//'     eigvals <- D_inv
//'   } else {
//'     eigvals <- D_noinv
//'   }
//'   return(rev_eig(eigvals, ED$vectors))
//' }
//'
//' # As armaRidgePAnyTarget
//' rRidgePAnyTarget <- function(S, tgt, l, invert = 2) {
//'   ED <- eigen(S - l*tgt, symmetric = TRUE)
//'   eigvals <- 0.5*ED$values
//'   sqroot <- sqrt(l + eigvals^2)
//'
//'   if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
//'     return(tgt)
//'   }
//'
//'   D_inv <- 1.0/(sqroot + eigvals)
//'   D_noinv <- (sqroot - eigvals)/l
//'
//'   if (invert == 2) {   # Determine to invert or not
//'     if (l > 1) {  # Generally, use inversion for "small" lambda
//'       invert = 0;
//'     } else {
//'       invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
//'     }
//'   }
//'
//'   if (invert == 1) {
//'     eigvals <- D_inv
//'   } else {
//'     eigvals <- D_noinv
//'   }
//'   return(rev_eig(eigvals, ED$vectors))
//' }
//'
//' rRidgeP <- function(S, tgt, l, invert = 2) {
//'   if (l == Inf) {
//'     return(tgt)
//'   }
//'   a <- tgt[1,1]
//'   if (tgt == diag(a, nrow(tgt))) {
//'     rRidgePScalarTarget(S, tgt, l, invert)
//'   } else {
//'     rRidgePAnyTarget(S, tgt, l, invert)
//'   }
//'
//' }
//'
//' # Contrasted to the straight forward implementations:
//' sqrtm <- function(X) { # Matrix square root
//'   ed <- eigen(X, symmetric = TRUE)
//'   rev_eig(sqrt(ed$values), ed$vectors)
//' }
//'
//' # Straight forward (Lemma 1)
//' ridgeP1 <- function(S, tgt, l) {
//'   solve(sqrtm( l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt) ) + 0.5*(S - l*tgt))
//' }
//'
//' # Straight forward  (Lemma 1 + remark 6 + 7)
//' ridgeP2 <- function(S, tgt, l) {
//'   1.0/l*(sqrtm(l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt)) - 0.5*(S - l*tgt))
//' }
//'
//' set.seed(1)
//' n <- 3
//' p <- 6
//' S <- covML(matrix(rnorm(p*n), n, p))
//' a <- 2.2
//' tgt <- diag(a, p)
//' l <- 1.21
//'
//' (A <- ridgeP1(S, tgt, l))
//' (B <- ridgeP2(S, tgt, l))
//'
//' (C  <- rags2ridges:::.armaRidgePAnyTarget(S, tgt, l))
//' (CR <-                   rRidgePAnyTarget(S, tgt, l))
//' (D  <- rags2ridges:::.armaRidgePScalarTarget(S, a, l))
//' (DR <-                   rRidgePScalarTarget(S, a, l))
//' (E  <- rags2ridges:::.armaRidgeP(S, tgt, l))
//'
//' # Check
//' equal <- function(x, y) {isTRUE(all.equal(x, y))}
//' stopifnot(equal(A, B) & equal(A, C) & equal(A, D) & equal(A, E))
//' stopifnot(equal(C, CR) & equal(D, DR))
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





//' Multivariate Gaussian simulation
//'
//' Fast simulation from multivariate Gaussian probability distribution.
//'
//' The \code{rmvnormal} function is copied from the \code{GMCM}-package. It is
//' similar to \code{rmvnorm} from the \code{mvtnorm}-package.
//'
//' @param n An \code{integer} giving the number of observations to be
//' simulated.
//' @param mu A \code{numeric} vector of dimension \eqn{p} giving the means of
//' normal distribution.
//' @param sigma A variance-covariance \code{matrix} of dimension \eqn{p} times
//' \eqn{p}.
//' @return Returns a \eqn{n} by \eqn{p} matrix of observations from a
//' multivariate normal distribution with the given mean \code{mu} and
//' covariance
//' @author Anders Ellern Bilgrau
//' @examples
//'
//' rmvnormal(n = 10, mu = 1:4, sigma = diag(4))
//'
//' @export rmvnormal
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
  const arma::mat L = arma::chol(sigma, "lower");
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
  const arma::mat L = arma::chol(inv(psi), "lower");
  arma::cube ans(p, p, n);
  for (int g = 0; g < n; ++g) {
    ans.slice(g) = inv(armaRWishartSingle(L, nu, p));
  }
  return ans;
}




////////////////////////////////////////////////////////////////////////////////
/* -----------------------------------------------------------------------------
SUPPORT FUNCTIONS FOR RIDGE PRECISION ESTIMATION WITH CHORDAL SUPPORT
INTENDED TO BE INCLUDED IN RAGS2RIDGES-PACKAGE WHEN DIRECT IMPORTATION FROM THAT
 SOURCE IS POSSIBLE
----------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(.armaRidgePchordalInit)]]
arma::mat armaRidgePchordalInitWorkhorse(arma::mat S, const double lambda, arma::mat target, std::string type, Rcpp::List Cliques, Rcpp::List Separators){
  ////////////////////////////////////////////////////////////////////////////
  // construct (guess of) ridge MLE of precision
  // various types (Alt, ArchI, ArchII) of the ridge estimor are implemented
  ////////////////////////////////////////////////////////////////////////////

  // start with null matrix
  arma::mat P = arma::zeros<arma::mat>(S.n_rows, S.n_rows);

  // first archetypal ridge estimator
  if (type == "ArchI"){
    arma::mat targetInv = arma::inv_sympd(target);
    // account for clique contributions
    for (int k = 0; k < Cliques.size(); ++k){
      arma::ivec slh = Cliques[k];
      slh = slh - 1;
      arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
      P.submat(cliqNodes, cliqNodes) =  P.submat(cliqNodes, cliqNodes) +
        arma::inv_sympd((1-lambda) * S.submat(cliqNodes, cliqNodes) + lambda * targetInv.submat(cliqNodes, cliqNodes));
    }
    // account for separator contributions
    for (int k = 0; k < Separators.size(); ++k){
      arma::ivec slh = Separators[k];
      if (slh.n_elem > 0){
        slh = slh - 1;
        arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
        P.submat(sepNodes, sepNodes) =  arma::symmatl(P.submat(sepNodes, sepNodes)) -
          arma::inv_sympd((1-lambda) * S.submat(sepNodes, sepNodes) + lambda * targetInv.submat(sepNodes, sepNodes));
      }
    }
  }

  // second archetypal ridge estimator
  if (type == "ArchII"){
    // account for clique contributions
    for (int k = 0; k < Cliques.size(); ++k){
      arma::ivec slh = Cliques[k];
      slh = slh - 1;
      arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
      P.submat(cliqNodes, cliqNodes) = arma::symmatl(P.submat(cliqNodes, cliqNodes)) +
        arma::inv_sympd(S.submat(cliqNodes, cliqNodes) + lambda * arma::eye(cliqNodes.n_elem, cliqNodes.n_elem));
    }
    // account for separator contributions
    for (int k = 0; k < Separators.size(); ++k){
      arma::ivec slh = Separators[k];
      if (slh.n_elem > 0){
        slh = slh - 1;
        arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
        P.submat(sepNodes, sepNodes) = arma::symmatl(P.submat(sepNodes, sepNodes)) -
          arma::inv_sympd(S.submat(sepNodes, sepNodes) + lambda * arma::eye(sepNodes.n_elem, sepNodes.n_elem));
      }
    }
  }

  // alternative ridge estimator
  if (type == "Alt"){
    // account for clique contributions
    for (int k = 0; k < Cliques.size(); ++k){
      arma::ivec slh = Cliques[k];
      slh = slh - 1;
      arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
      if (cliqNodes.n_elem > 1){ P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) +
        armaRidgeP(S.submat(cliqNodes, cliqNodes), target.submat(cliqNodes, cliqNodes), lambda);
      }
      if (cliqNodes.n_elem == 1){
        arma::mat slh2 = S.submat(cliqNodes, cliqNodes) - lambda * target.submat(cliqNodes, cliqNodes);
        double slh3 = slh2(0,0);
        P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) + 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
      }
    }
    // account for separator contributions
    for (int k = 0; k < Separators.size(); ++k){
      arma::ivec slh = Separators[k];
      if (slh.n_elem > 0){
        slh = slh - 1;
        arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
        if (sepNodes.n_elem > 1){ P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) -
          armaRidgeP(S.submat(sepNodes, sepNodes), target.submat(sepNodes, sepNodes), lambda);
        }
        if (sepNodes.n_elem == 1){
          arma::mat slh2 = S.submat(sepNodes, sepNodes) - lambda * target.submat(sepNodes, sepNodes);
          double slh3 = slh2(0,0);
          P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) - 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
        }
      }
    }
  }

  // return initial ridge guess for chordal precision matrix
  return(P);
}

// [[Rcpp::export(.armaPenLLreparPforNLM)]]
Rcpp::NumericVector armaPenLLreparPforNLM(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
  /////////////////////////////////////////////////////////////////////////////////////////
  // ridge penalized log-likelihood (and its gradient) for the reparametrized precision matrix.
  // for the reparametrization refer Dahl et al. (2005).
  // passed on to the 'nlm' optimization function
  /////////////////////////////////////////////////////////////////////////////////////////

  // construct precision matrix from alternative parametrization
  const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

  // return (minus) penalized log-likelihood
  Rcpp::NumericVector penLL(1);
  penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));

  // return gradient of (minus) penalized log-likelihood: nonzero elements only
  int p = S.n_rows;
  const arma::mat gradMat = 2 * (S - arma::inv_sympd(P) + lambda * (P - target));
  arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
  arma::vec grad = gradMat.elem(elemID);
  penLL.attr("gradient") = grad;

  return penLL;
}

// [[Rcpp::export(.armaPenLLreparP)]]
const double armaPenLLreparP(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
  /////////////////////////////////////////////////////////////////////////////////////////
  // ridge penalized log-likelihood for the reparametrized precision matrix.
  // for the reparametrization refer Dahl et al. (2005)
  // passed on to the 'optim' optimization function
  /////////////////////////////////////////////////////////////////////////////////////////

  // construct precision matrix from alternative parametrization
  const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

  // return (minus) penalized log-likelihood
  // double logDetP; double detPsign;
  // arma::log_det(logDetP, detPsign, P);
  // const double penLL = -logDetP + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));
  const double penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));
  return penLL;
}

// [[Rcpp::export(.armaPenLLreparPgrad)]]
arma::vec armaPenLLreparPgrad(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
  /////////////////////////////////////////////////////////////////////////////////////////
  // gradient of ridge penalized log-likelihood for the reparametrized precision matrix.
  // for the reparametrization refer Dahl et al. (2005)
  // passed on to the 'optim' optimization function
  /////////////////////////////////////////////////////////////////////////////////////////

  // construct precision matrix from alternative parametrization
  const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

  // return gradient of (minus) penalized log-likelihood: nonzero elements only
  int p = S.n_rows;
  const arma::mat gradMat = 2 * (arma::symmatl(S) - arma::symmatl(arma::inv_sympd(P)) + lambda * (arma::symmatl(P) - arma::symmatl(target)));
  arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
  return gradMat.elem(elemID);
}

// [[Rcpp::export(.armaPenLLreparGradArchI)]]
arma::vec armaPenLLreparGradArchI(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
  /////////////////////////////////////////////////////////////////////////////////////////
  // gradient of 'ArchI'-ridge penalized log-likelihood for the reparametrized precision matrix.
  // for the reparametrization refer Dahl et al. (2005)
  // passed on to the 'optim' optimization function
  /////////////////////////////////////////////////////////////////////////////////////////

  // construct precision matrix from alternative parametrization
  const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

  // return gradient of (minus) penalized log-likelihood: nonzero elements only
  int p = S.n_rows;
  const arma::mat gradMat = 2 * ((1-lambda) * arma::symmatl(S) + lambda * arma::symmatl(arma::inv_sympd(target)) - arma::symmatl(arma::inv_sympd(P)));
  arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
  return gradMat.elem(elemID);
}

// [[Rcpp::export(.armaPenLLreparGradArchII)]]
arma::vec armaPenLLreparGradArchII(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
  /////////////////////////////////////////////////////////////////////////////////////////
  // gradient of 'ArchII-ridge penalized log-likelihood for the reparametrized precision matrix.
  // for the reparametrization refer Dahl et al. (2005)
  // passed on to the 'optim' optimization function
  /////////////////////////////////////////////////////////////////////////////////////////

  // construct precision matrix from alternative parametrization
  const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

  // return gradient of (minus) penalized log-likelihood: nonzero elements only
  int p = S.n_rows;
  const arma::mat gradMat = 2 * (arma::symmatl(S) + lambda * arma::eye(p,p) - arma::inv_sympd(P));
  arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
  return gradMat.elem(elemID);
}


