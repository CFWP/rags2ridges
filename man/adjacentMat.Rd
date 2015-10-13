\name{adjacentMat}
\alias{adjacentMat}
\title{
Transform real matrix into an adjacency matrix
}
\description{
Function that transforms a real matrix into an adjacency matrix. Intended use: Turn sparsified precision matrix into an adjacency matrix for undirected graphical representation.
}
\usage{
adjacentMat(M, diag = FALSE)
}
\arguments{
  \item{M}{
(Possibly sparsified precision) \code{matrix}.
}
  \item{diag}{
A \code{logical} indicating if the diagonal elements should be retained.
}
}
\value{
Function returns an adjacency \code{matrix}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeP}}, \code{\link{covML}}, \code{\link{sparsify}}, \code{\link{edgeHeat}}, \code{\link{Ugraph}}
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

## Obtain sparsified partial correlation matrix
PC0 <- sparsify(P, threshold = "localFDR", FDRcut = .8)

## Obtain adjacency matrix
adjacentMat(PC0$sparsePrecision)
}
