\name{GGMpathStats.fused}
\alias{GGMpathStats.fused}
\title{
Fused gaussian graphical model node pair path statistics
}
\description{
A simple wrapper for \code{\link{GGMpathStats}}.
}
\usage{
GGMpathStats.fused(sparsePList, ...) 
}
\arguments{
  \item{sparsePList}{
    A \code{list} of sparsified precision matrices.
}
  \item{\dots}{
    Arguments passed to \code{\link{GGMpathStats}}.
}
}
\value{
A \code{list} of path stats.
}
\note{
The function currently fails if no paths are preset in one of the groups.
}
\author{
Anders Ellern Bilgrau
}
\seealso{
\code{\link{GGMpathStats}}
}
\examples{
## Obtain some (high-dimensional) data
set.seed(10)
ns <- c(10, 11)
SList <- createS(ns, 7)
TList <- default.target.fused(SList, ns)

## Obtain regularized precision and sparsify
PList <- ridgeS.fused(SList, ns, TList, lambda = 0.3, lambdaF = 1.5)
sparsePList <- sparsify.fused(PList, threshold = "absValue", absValueCut = 0.20)
sparsePList <- lapply(sparsePList, "[[", "sparsePrecision")

## Obtain information on mediating and moderating paths between nodes 14 and 23
res <- GGMpathStats.fused(sparsePList, node1 = 1, node2 = 2, graph = FALSE)
}
