\name{pruneMatrix}
\alias{pruneMatrix}
\title{
Prune square matrix to those variables having nonzero entries
}
\description{
Convenience function that prunes a square matrix to those variables (features) having nonzero row (column) entries (i.e., to features implied in graphical connections).
}
\usage{
pruneMatrix(M)
}
\arguments{
  \item{M}{
(Possibly sparsified) square \code{matrix}.
}
}
\value{
A pruned \code{matrix}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)

## Determine support regularized standardized precision under optimal penalty
PC0 <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor

## Prune sparsified partial correlation matrix
PC0P <- pruneMatrix(PC0)
}
