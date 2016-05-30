\name{.armaRidgeP}
\alias{.armaRidgeP}
\alias{.armaRidgePAnyTarget}
\alias{.armaRidgePScalarTarget}
\alias{armaRidgeP}
\alias{armaRidgePAnyTarget}
\alias{armaRidgePScalarTarget}
\title{
Core ridge precision estimators
}
\description{
This is the \R interface to the \code{C++} implementations of the ridge precision estimators.
They perform core computations for many other functions.
}
\usage{
# General usage
.armaRidgeP(S, target, lambda, invert = 2L)

# Arbitrary target estimator
.armaRidgePAnyTarget(S, target, lambda, invert = 2L)

# Scalar target (alpha*I) estimator
.armaRidgePScalarTarget(S, alpha, lambda, invert = 2L)
}
\arguments{
  \item{S}{
  A sample covariance \code{matrix}.
  }
  \item{target}{
  A \code{numeric} symmetric positive definite target \code{matrix} of the same
  size as \code{S}.
  }
  \item{alpha}{
  The diagonal value on the scalar target matrix. A \code{double} of length 1.
  }
  \item{lambda}{
  The ridge penalty. A \code{double} of length 1.
  }
  \item{invert}{
  An \code{integer}. Should the estimate be computed using inversion?
  Permitted values are \code{0L}, \code{1L}, and \code{2L} meaning "no", "yes",
  and "automatic" (default).
  }
}
\details{
  The functions are R-interfaces to low-level \code{C++} implementations
  of the ridge estimators in the reference below
  (cf. Lemma 1, Remark 6, Remark 7, and section 4.1 therein).

  \code{.armaRidgeP} is simply a wrapper (on the C++ side) for
  \code{.armaRidgePAnyTarget} and \code{.armaRidgePScalarTarget} which are
  the estimators for arbitrary and scalar targets, respectively.
  The \code{invert} argument of the functions indicates if the computation
  uses matrix inversion or not.

  Essentially, the functions compute
  \deqn{
    \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
    \left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
    \lambda_{a}\mathbf{T})^{2}\right]^{1/2} + \frac{1}{2}(\mathbf{S} -
    \lambda_{a}\mathbf{T})\right\}^{-1},
  }{%
    \Omega(\lambda) =
    \{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
    + 1/2(S - \lambda T)\}^{-1},
  }
  if \code{invert == 1} or
  \deqn{
    \hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
    \frac{1}{\lambda}\left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
    \lambda_{a}\mathbf{T})^{2}\right]^{1/2} - \frac{1}{2}(\mathbf{S} -
    \lambda_{a}\mathbf{T})\right\}
  }{%
    \Omega(\lambda) =
    1/{[\lambda I + 1/4(S - \lambda T)^2 ]^{1/2}
    - 1/2(S - \lambda T)},
  }
  if \code{invert == 0} using appropriate eigenvalue decompositions.
  See the \R implementations in the example section below.
}
\value{
  Returns a symmetric positive definite \code{matrix} of the same size as
  \code{S}.
}
\references{
  van Wieringen, W.N. & Peeters, C.F.W. (To appear).
  Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
  Computational Statistics & Data Analysis.
  Available as arXiv:1403.0904v3 [stat.ME].
}
\author{
  Anders Ellern Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\section{Warning}{
  The functions themselves perform no checks on the input. Correct input should
  be ensured by wrappers.
}
\seealso{
  Used as a backbone in \code{\link{ridgeP}}, \code{\link{ridgeP.fused}}, etc.
}
\examples{
S <- createS(n = 3, p = 4)
tgt <- diag(4)
rags2ridges:::.armaRidgeP(S, tgt, 1.2)

rags2ridges:::.armaRidgePAnyTarget(S, tgt, 1.2)
rags2ridges:::.armaRidgePScalarTarget(S, 1, 1.2)


################################################################################
# The C++ estimators essentially amount to the following functions.
################################################################################

rev_eig <- function(evalues, evectors) { # "Reverse" eigen decomposition
  evectors \%*\% diag(evalues) \%*\% t(evectors)
}

# R implementations

# As armaRidgePScalarTarget Inv
rRidgePScalarTarget <- function(S, a, l, invert = 2) {
  ED <- eigen(S, symmetric = TRUE)
  eigvals <- 0.5*(ED$values - l*a)
  sqroot <- sqrt(l + eigvals^2)

  if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
    return(diag(a, nrow(S)))
  }

  D_inv <- 1.0/(sqroot + eigvals)
  D_noinv <- (sqroot - eigvals)/l

  if (invert == 2) {   # Determine to invert or not
    if (l > 1) {  # Generally, use inversion for "small" lambda
      invert = 0;
    } else {
      invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
    }
  }

  if (invert) {
    eigvals <- D_inv
  } else {
    eigvals <- D_noinv
  }
  return(rev_eig(eigvals, ED$vectors))
}

# As armaRidgePAnyTarget
rRidgePAnyTarget <- function(S, tgt, l, invert = 2) {
  ED <- eigen(S - l*tgt, symmetric = TRUE)
  eigvals <- 0.5*ED$values
  sqroot <- sqrt(l + eigvals^2)

  if (l > 1e6 && (any(!is.finite(eigvals)) || any(!is.finite(sqroot)))) {
    return(tgt)
  }

  D_inv <- 1.0/(sqroot + eigvals)
  D_noinv <- (sqroot - eigvals)/l

  if (invert == 2) {   # Determine to invert or not
    if (l > 1) {  # Generally, use inversion for "small" lambda
      invert = 0;
    } else {
      invert <- ifelse(any(!is.finite(D_inv)), 0, 1)
    }
  }

  if (invert == 1) {
    eigvals <- D_inv
  } else {
    eigvals <- D_noinv
  }
  return(rev_eig(eigvals, ED$vectors))
}

rRidgeP <- function(S, tgt, l, invert = 2) {
  if (l == Inf) {
    return(tgt)
  }
  a <- tgt[1,1]
  if (tgt == diag(a, nrow(tgt))) {
    rRidgePScalarTarget(S, tgt, l, invert)
  } else {
    rRidgePAnyTarget(S, tgt, l, invert)
  }

}

# Contrasted to the straight forward implementations:
sqrtm <- function(X) { # Matrix square root
  ed <- eigen(X, symmetric = TRUE)
  rev_eig(sqrt(ed$values), ed$vectors)
}

# Straight forward (Lemma 1)
ridgeP1 <- function(S, tgt, l) {
  solve(sqrtm( l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt) ) + 0.5*(S - l*tgt))
}

# Straight forward  (Lemma 1 + remark 6 + 7)
ridgeP2 <- function(S, tgt, l) {
  1.0/l*(sqrtm(l*diag(nrow(S)) + 0.25*crossprod(S - l*tgt)) - 0.5*(S - l*tgt))
}

set.seed(1)
n <- 3
p <- 6
S <- covML(matrix(rnorm(p*n), n, p))
a <- 2.2
tgt <- diag(a, p)
l <- 1.21

(A <- ridgeP1(S, tgt, l))
(B <- ridgeP2(S, tgt, l))

(C  <- rags2ridges:::.armaRidgePAnyTarget(S, tgt, l))
(CR <-                   rRidgePAnyTarget(S, tgt, l))
(D  <- rags2ridges:::.armaRidgePScalarTarget(S, a, l))
(DR <-                   rRidgePScalarTarget(S, a, l))
(E  <- rags2ridges:::.armaRidgeP(S, tgt, l))

# Check
equal <- function(x, y) {isTRUE(all.equal(x, y))}
stopifnot(equal(A, B) & equal(A, C) & equal(A, D) & equal(A, E))
stopifnot(equal(C, CR) & equal(D, DR))

}
\keyword{internal}
