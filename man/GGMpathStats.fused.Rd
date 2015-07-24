\name{GGMpathStats.fused}
\alias{GGMpathStats.fused}
\title{
Fused gaussian graphical model node pair path statistics
}
\description{
A simple wrapper for \code{\link{GGMpathStats}}.
}
\usage{
GGMpathStats.fused(sparsePlist, ...)
}
\arguments{
  \item{sparsePlist}{
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
Slist <- createS(ns, p = 7, topology = "banded")
Tlist <- default.target.fused(Slist, ns)

## Obtain regularized precision and sparsify
Plist <- ridgeP.fused(Slist, ns, Tlist, lambda = c(1, 1.6))
sparsePlist <- sparsify.fused(Plist, threshold = "absValue", absValueCut = 0.20)
SPlist <- lapply(sparsePlist, "[[", "sparsePrecision")

## Obtain information on mediating and moderating paths between nodes 14 and 23
res <- GGMpathStats.fused(SPlist, node1 = 3, node2 = 4, graph = FALSE)

}
