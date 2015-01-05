\name{symm}
\alias{symm}
\title{
Symmetrize matrix
}
\description{
Function that symmetrizes matrices.
}
\usage{
symm(M)
}
\arguments{
  \item{M}{
(In numeric ideality symmetric) square \code{matrix}.
}
}
\details{
Large objects that are symmetric sometimes fail to be recognized as such by R due to
rounding under machine precision. This function symmetrizes for computational purposes
matrices that are symmetric in numeric ideality.
}
\value{
A symmetric \code{matrix}. 
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT <- optPenalty.LOOCV(X, 10, 30, 10, target = diag(diag(1/covML(X))))

## Check symmetry
## OPT$optPrec is symmetric by definition
## But is not recognized as such due to rounding peculiarities
isSymmetric(OPT$optPrec)

## Symmetrize
symm(OPT$optPrec)
}