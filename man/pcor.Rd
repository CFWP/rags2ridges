\name{pcor}
\alias{pcor}
\title{
Compute partial correlation matrix or standardized precision matrix
}
\description{
Function computing the partial correlation matrix or standardized precision matrix from an input precision matrix.
}
\usage{
pcor(P, pc = TRUE)
}
\arguments{
  \item{P}{
(Possibly regularized) precision \code{matrix}.
}
  \item{pc}{
A \code{logical} indicating if the partial correlation matrix should be computed.
}
}
\details{
The function assumes that the input \code{matrix} is a precision matrix.
If \code{pc = FALSE} the standardized precision matrix, rather than the partial correlation matrix, is given as the output value.
The standardized precision matrix is equal to the partial correlation matrix up to the sign of off-diagonal entries.
}
\value{
A partial correlation \code{matrix} or a standardized precision \code{matrix}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeP}}, \code{\link{covML}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Obtain regularized precision matrix
P <- ridgeP(Cx, lambda = 10, type = "Alt")

## Obtain partial correlation matrix
pcor(P)
}
