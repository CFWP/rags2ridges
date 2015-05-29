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
This is the R interface to the \code{C++} implementations of the
ridge precision estimators. They serves as core functions for many other
functions.
}
\usage{
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
  A \code{numeric} symmetric positive definite target \code{matrix} the same
  size as \code{S}.
  }
  \item{alpha}{
  The diagonal value on the scalar target matrix. A \code{double} of length 1.
  }
  \item{lambda}{
  The ridge penalty. A \code{double} of length 1.
  }
  \item{invert}{
  An \code{integer}. Should the estimate be compute using inversion?
  Permitted values are 0 = "no", 1 = "yes", and 2 = "automatic" (default).
  }
}
\details{
  The function are interface the low-level \code{C++} implementations
  of the ridge estimators in the reference below.
  Cf. Lemma 1, Remark 6, Remark 7, and section 4.1 herein.

  \code{.armaRidgeP} is simply a wrapper (on the C++ side) for
  \code{.armaRidgePAnyTarget} and \code{.armaRidgePScalarTarget} which are
  the estimators for arbitrary and scalar targets, respectively.
  The \code{invert} argument of the functions indicate that the computation
  uses matrix inversion or not.
}
\value{
  Returns a symmetric positive definite \code{matrix} of the same size as \code{S}.
}
\references{
  van Wieringen, W.N. and Peeters, C.F.W. (2014). Ridge Estimation of Inverse
  Covariance Matrices from High-Dimensional Data. arXiv:1403.0904 [stat.ME].
}
\author{
  Anders Ellern Bilgrau
  Carel Peeters
  Wessel van Wieringen
}
\section{Warning}{
  The functions themselves perform very few checks on the input.
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





#
# The estimators essentially amount C++ versions of the following functions.
#

rev_eig <- function(evalues, evectors) { # "Reverse" eigen decomposition
  evectors \%*\% diag(evalues) \%*\% t(evectors)
}

# R implementations

# As armaRidgePScalarTarget Inv
ridgePScalarTarget <- function(S, a, l, invert = 2) {
  ED <- eigen(S)
  tmp <- ED$values - l*a
  d <- (sqrt(l + 0.25*(tmp)^2) - 0.5*tmp)

  if (invert == 2)
    invert <- ifelse(all(d == 0), 1, 0)

  if (invert) {
    eigvals <- 1.0/(sqrt(l + 0.25*(tmp)^2) + 0.5*tmp)
  } else {
    eigvals <- (1.0/l)*d
  }
  return(rev_eig(eigvals, ED$vectors))
}

# As armaRidgePAnyTarget
ridgePAnyTarget <- function(S, tgt, l, invert = 2) {
  ED <- eigen(S - l*tgt)

  d <- (sqrt(l + 0.25*(ED$values)^2) - 0.5*ED$values)

  if (invert == 2)
    invert <- ifelse(all(d == 0), 1, 0)

  if (invert == 1) {
    eigvals <- 1.0/(sqrt(l + 0.25*(ED$values^2)) + 0.5*(ED$values))
  } else {
    eigvals <- (1.0/l)*d
  }
  return(rev_eig(eigvals, ED$vectors))
}

# Contrasted to the straight forward implementations:
sqrtm <- function(X) { # Matrix square root
  ed <- eigen(X)
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

(AA <- ridgeP1(S, tgt, l))
(BB <- ridgeP2(S, tgt, l))

(CC <- rags2ridges:::.armaRidgePAnyTarget(S, tgt, l))
(CCR<- ridgePAnyTarget(S, tgt, l))
(DD <- rags2ridges:::.armaRidgePScalarTarget(S, a, l))
(DDR<- ridgePScalarTarget(S, a, l))
(EE <- rags2ridges:::.armaRidgeP(S, tgt, l))

# Check
equal <- function(x, y) {isTRUE(all.equal(x, y))}
equal(AA, BB) & equal(AA, CC) & equal(AA, DD) & equal(AA, EE)
equal(CC, CCR) & equal(DD, DDR)

}
\keyword{internal}
