# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Evaluate the (penalized) (fused) likelihood
#'
#' Functions that evaluate the (penalized) (fused) likelihood.
#'
#' @aliases NLL PNLL NLL.fused PNLL.fused
#' @param S,Slist A (list of) positive semi definite sample covariance
#' matrices.
#' @param P,Plist A (list of) positive definite precision matrices.
#' @param T,Tlist A (list of) positive definite target matrices.
#' @param ns A \code{numeric} of sample sizes.
#' @param lambda A \code{numeric} penalty parameter. For the \code{.fused}
#' functions, this is a penalty \code{matrix}.
#' @return A single number.
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#' Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{ridgeP.fused}}
#' @examples
#'
#' ns <- c(4,5)
#' Slist <- createS(n = ns, p = 5)
#' Plist <- list(diag(5), diag(2,5))
#' Tlist <- list(diag(5), diag(5))
#'
#' NLL(Slist[[1]], Plist[[1]])
#' PNLL(Slist[[1]], Plist[[1]], Tlist[[1]], lambda = 1)
#' NLL.fused(Slist, Plist, ns)
#' PNLL.fused(Slist, Plist, ns, Tlist, lambda = diag(2))
#'
#' @export
NLL <- function(S, P) {
    .Call('_rags2ridges_NLL', PACKAGE = 'rags2ridges', S, P)
}

#' @rdname NLL
#' @export
PNLL <- function(S, P, T, lambda) {
    .Call('_rags2ridges_PNLL', PACKAGE = 'rags2ridges', S, P, T, lambda)
}

#' @rdname NLL
#' @export
NLL.fused <- function(Slist, Plist, ns) {
    .Call('_rags2ridges_NLL_fused', PACKAGE = 'rags2ridges', Slist, Plist, ns)
}

#' @rdname NLL
#' @export
PNLL.fused <- function(Slist, Plist, ns, Tlist, lambda) {
    .Call('_rags2ridges_PNLL_fused', PACKAGE = 'rags2ridges', Slist, Plist, ns, Tlist, lambda)
}

.armaPooledS <- function(Slist, ns, mle = 0L) {
    .Call('_rags2ridges_armaPooledS', PACKAGE = 'rags2ridges', Slist, ns, mle)
}

.armaPooledP <- function(Plist, ns, mle = 0L) {
    .Call('_rags2ridges_armaPooledP', PACKAGE = 'rags2ridges', Plist, ns, mle)
}

.armaEigShrink <- function(dVec, lambda, cons = 0) {
    .Call('_rags2ridges_armaEigShrink', PACKAGE = 'rags2ridges', dVec, lambda, cons)
}

.armaEigShrinkAnyTarget <- function(S, target, lambda) {
    .Call('_rags2ridges_armaEigShrinkAnyTarget', PACKAGE = 'rags2ridges', S, target, lambda)
}

.armaEigShrinkArchI <- function(dVec, lambda, cons) {
    .Call('_rags2ridges_armaEigShrinkArchI', PACKAGE = 'rags2ridges', dVec, lambda, cons)
}

.armaRidgePAnyTarget <- function(S, target, lambda, invert = 2L) {
    .Call('_rags2ridges_armaRidgePAnyTarget', PACKAGE = 'rags2ridges', S, target, lambda, invert)
}

.armaRidgePScalarTarget <- function(S, alpha, lambda, invert = 2L) {
    .Call('_rags2ridges_armaRidgePScalarTarget', PACKAGE = 'rags2ridges', S, alpha, lambda, invert)
}

#' Core ridge precision estimators
#'
#' This is the interface to the \code{C++} implementations of the ridge
#' precision estimators. They perform core computations for many other
#' functions.
#'
#' @details
#' The functions are R-interfaces to low-level \code{C++} implementations
#' of the ridge estimators in the reference below
#' (cf. Lemma 1, Remark 6, Remark 7, and section 4.1 therein).
#'
#' \code{.armaRidgeP} is simply a wrapper (on the C++ side) for
#' \code{.armaRidgePAnyTarget} and \code{.armaRidgePScalarTarget} which are
#' the estimators for arbitrary and scalar targets, respectively.
#' The \code{invert} argument of the functions indicates if the computation
#' uses matrix inversion or not.
#'
#' Essentially, the functions compute
#'     \deqn{
#'       \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
#'         \left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
#'         \lambda_{a}\mathbf{T})^{2}\right]^{1/2} + \frac{1}{2}(\mathbf{S} -
#'         \lambda_{a}\mathbf{T})\right\}^{-1},
#'     }{%
#'       \Omega(\lambda) =
#'       \{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
#'         + 1/2(S - \lambda T)\}^{-1},
#'     }
#' if \code{invert == 1} or
#'     \deqn{
#'       \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
#'         \frac{1}{\lambda}\left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
#'                  \lambda_{a}\mathbf{T})^{2}\right]^{1/2} - \frac{1}{2}(\mathbf{S} -
#'         \lambda_{a}\mathbf{T})\right\}
#'     }{%
#'       \Omega(\lambda) =
#'       1/{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
#'         - 1/2(S - \lambda T)},
#'     }
#' if \code{invert == 0} using appropriate eigenvalue decompositions.
#' See the \R implementations in the example section below.
#'
#' @aliases .armaRidgeP .armaRidgePAnyTarget .armaRidgePScalarTarget armaRidgeP
#' armaRidgePAnyTarget armaRidgePScalarTarget
#' @param S A sample covariance \code{matrix}.
#' @param target A \code{numeric} symmetric positive definite target
#' \code{matrix} of the same size as \code{S}.
#' @param lambda The ridge penalty. A \code{double} of length 1.
#' @param invert An \code{integer}. Should the estimate be computed using
#' inversion?  Permitted values are \code{0L}, \code{1L}, and \code{2L} meaning
#' "no", "yes", and "automatic" (default).
#' @return Returns a symmetric positive definite \code{matrix} of the same size
#' as \code{S}.
#' @section Warning: The functions themselves perform no checks on the input.
#' Correct input should be ensured by wrappers.
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso Used as a backbone in \code{\link{ridgeP}},
#' \code{\link{ridgeP.fused}}, etc.
#' @references van Wieringen, W.N. & Peeters, C.F.W. (2016).  Ridge Estimation
#' of Inverse Covariance Matrices from High-Dimensional Data, Computational
#' Statistics & Data Analysis, vol. 103: 284-303.  Also available as
#' arXiv:1403.0904v3 [stat.ME].
#' @keywords internal
#' @examples
#'
#' S <- createS(n = 3, p = 4)
#' tgt <- diag(4)
#' rags2ridges:::.armaRidgeP(S, tgt, 1.2)
#'
#' rags2ridges:::.armaRidgePAnyTarget(S, tgt, 1.2)
#' rags2ridges:::.armaRidgePScalarTarget(S, 1, 1.2)
#'
#'
#' ################################################################################
#' # The C++ estimators essentially amount to the following functions.
#' ################################################################################
#'
#' rev_eig <- function(evalues, evectors) { # "Reverse" eigen decomposition
#'   evectors %*% diag(evalues) %*% t(evectors)
#' }
#'
#' # R implementations
#'
#' # As armaRidgePScalarTarget Inv
#' rRidgePScalarTarget <- function(S, a, l, invert = 2) {
#'   ED <- eigen(S, symmetric = TRUE)
#'   eigvals <- 0.5*(ED$values - l*a)
#'   sqroot <- sqrt(l + eigvals^2)
#'
#'   if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
#'     return(diag(a, nrow(S)))
#'   }
#'
#'   D_inv <- 1.0/(sqroot + eigvals)
#'   D_noinv <- (sqroot - eigvals)/l
#'
#'   if (invert == 2) {   # Determine to invert or not
#'     if (l > 1) {  # Generally, use inversion for "small" lambda
#'       invert = 0;
#'     } else {
#'       invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
#'     }
#'   }
#'
#'   if (invert) {
#'     eigvals <- D_inv
#'   } else {
#'     eigvals <- D_noinv
#'   }
#'   return(rev_eig(eigvals, ED$vectors))
#' }
#'
#' # As armaRidgePAnyTarget
#' rRidgePAnyTarget <- function(S, tgt, l, invert = 2) {
#'   ED <- eigen(S - l*tgt, symmetric = TRUE)
#'   eigvals <- 0.5*ED$values
#'   sqroot <- sqrt(l + eigvals^2)
#'
#'   if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
#'     return(tgt)
#'   }
#'
#'   D_inv <- 1.0/(sqroot + eigvals)
#'   D_noinv <- (sqroot - eigvals)/l
#'
#'   if (invert == 2) {   # Determine to invert or not
#'     if (l > 1) {  # Generally, use inversion for "small" lambda
#'       invert = 0;
#'     } else {
#'       invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
#'     }
#'   }
#'
#'   if (invert == 1) {
#'     eigvals <- D_inv
#'   } else {
#'     eigvals <- D_noinv
#'   }
#'   return(rev_eig(eigvals, ED$vectors))
#' }
#'
#' rRidgeP <- function(S, tgt, l, invert = 2) {
#'   if (l == Inf) {
#'     return(tgt)
#'   }
#'   a <- tgt[1,1]
#'   if (tgt == diag(a, nrow(tgt))) {
#'     rRidgePScalarTarget(S, tgt, l, invert)
#'   } else {
#'     rRidgePAnyTarget(S, tgt, l, invert)
#'   }
#'
#' }
#'
#' # Contrasted to the straight forward implementations:
#' sqrtm <- function(X) { # Matrix square root
#'   ed <- eigen(X, symmetric = TRUE)
#'   rev_eig(sqrt(ed$values), ed$vectors)
#' }
#'
#' # Straight forward (Lemma 1)
#' ridgeP1 <- function(S, tgt, l) {
#'   solve(sqrtm( l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt) ) + 0.5*(S - l*tgt))
#' }
#'
#' # Straight forward  (Lemma 1 + remark 6 + 7)
#' ridgeP2 <- function(S, tgt, l) {
#'   1.0/l*(sqrtm(l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt)) - 0.5*(S - l*tgt))
#' }
#'
#' set.seed(1)
#' n <- 3
#' p <- 6
#' S <- covML(matrix(rnorm(p*n), n, p))
#' a <- 2.2
#' tgt <- diag(a, p)
#' l <- 1.21
#'
#' (A <- ridgeP1(S, tgt, l))
#' (B <- ridgeP2(S, tgt, l))
#'
#' (C  <- rags2ridges:::.armaRidgePAnyTarget(S, tgt, l))
#' (CR <-                   rRidgePAnyTarget(S, tgt, l))
#' (D  <- rags2ridges:::.armaRidgePScalarTarget(S, a, l))
#' (DR <-                   rRidgePScalarTarget(S, a, l))
#' (E  <- rags2ridges:::.armaRidgeP(S, tgt, l))
#'
#' # Check
#' equal <- function(x, y) {isTRUE(all.equal(x, y))}
#' stopifnot(equal(A, B) & equal(A, C) & equal(A, D) & equal(A, E))
#' stopifnot(equal(C, CR) & equal(D, DR))
.armaRidgeP <- function(S, target, lambda, invert = 2L) {
    .Call('_rags2ridges_armaRidgeP', PACKAGE = 'rags2ridges', S, target, lambda, invert)
}

.armaFusedUpdateI <- function(g0, Plist, Slist, Tlist, ns, lambda) {
    .Call('_rags2ridges_armaFusedUpdateI', PACKAGE = 'rags2ridges', g0, Plist, Slist, Tlist, ns, lambda)
}

.armaFusedUpdateII <- function(g0, Plist, Slist, Tlist, ns, lambda) {
    .Call('_rags2ridges_armaFusedUpdateII', PACKAGE = 'rags2ridges', g0, Plist, Slist, Tlist, ns, lambda)
}

.armaFusedUpdateIII <- function(g0, Plist, Slist, Tlist, ns, lambda) {
    .Call('_rags2ridges_armaFusedUpdateIII', PACKAGE = 'rags2ridges', g0, Plist, Slist, Tlist, ns, lambda)
}

.armaRidgeP.fused <- function(Slist, ns, Tlist, lambda, Plist, maxit = 100L, eps = 1e-7, relative = TRUE, verbose = FALSE) {
    .Call('_rags2ridges_armaRidgeP_fused', PACKAGE = 'rags2ridges', Slist, ns, Tlist, lambda, Plist, maxit, eps, relative, verbose)
}

#' Multivariate Gaussian simulation
#'
#' Fast simulation from multivariate Gaussian probability distribution.
#'
#' The \code{rmvnormal} function is copied from the \code{GMCM}-package. It is
#' similar to \code{rmvnorm} from the \code{mvtnorm}-package.
#'
#' @param n An \code{integer} giving the number of observations to be
#' simulated.
#' @param mu A \code{numeric} vector of dimension \eqn{p} giving the means of
#' normal distribution.
#' @param sigma A variance-covariance \code{matrix} of dimension \eqn{p} times
#' \eqn{p}.
#' @return Returns a \eqn{n} by \eqn{p} matrix of observations from a
#' multivariate normal distribution with the given mean \code{mu} and
#' covariance
#' @author Anders Ellern Bilgrau
#' @examples
#'
#' rmvnormal(n = 10, mu = 1:4, sigma = diag(4))
#'
#' @export rmvnormal
rmvnormal <- function(n, mu, sigma) {
    .Call('_rags2ridges_rmvnormal', PACKAGE = 'rags2ridges', n, mu, sigma)
}

.armaRWishart <- function(n, sigma, nu) {
    .Call('_rags2ridges_armaRWishart', PACKAGE = 'rags2ridges', n, sigma, nu)
}

.armaRInvWishart <- function(n, psi, nu) {
    .Call('_rags2ridges_armaRInvWishart', PACKAGE = 'rags2ridges', n, psi, nu)
}

.armaRidgePchordalInit <- function(S, lambda, target, type, Cliques, Separators) {
    .Call('_rags2ridges_armaRidgePchordalInitWorkhorse', PACKAGE = 'rags2ridges', S, lambda, target, type, Cliques, Separators)
}

.armaPenLLreparPforNLM <- function(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC) {
    .Call('_rags2ridges_armaPenLLreparPforNLM', PACKAGE = 'rags2ridges', x, E1, E2, S, lambda, target, nonzerosR, nonzerosC)
}

.armaPenLLreparP <- function(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC) {
    .Call('_rags2ridges_armaPenLLreparP', PACKAGE = 'rags2ridges', x, E1, E2, S, lambda, target, nonzerosR, nonzerosC)
}

.armaPenLLreparPgrad <- function(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC) {
    .Call('_rags2ridges_armaPenLLreparPgrad', PACKAGE = 'rags2ridges', x, E1, E2, S, lambda, target, nonzerosR, nonzerosC)
}

.armaPenLLreparGradArchI <- function(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC) {
    .Call('_rags2ridges_armaPenLLreparGradArchI', PACKAGE = 'rags2ridges', x, E1, E2, S, lambda, target, nonzerosR, nonzerosC)
}

.armaPenLLreparGradArchII <- function(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC) {
    .Call('_rags2ridges_armaPenLLreparGradArchII', PACKAGE = 'rags2ridges', x, E1, E2, S, lambda, target, nonzerosR, nonzerosC)
}

