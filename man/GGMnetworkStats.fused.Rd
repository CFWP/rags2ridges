\name{GGMnetworkStats.fused}
\alias{GGMnetworkStats.fused}
\title{
Gaussian graphical model network statistics
}
\description{
Compute various network statistics from a \code{list} sparse precision matrices. The sparse precision matrix is taken to represent the conditional independence graph of a Gaussian graphical model.
This function is a simple wrapper for \code{\link{GGMnetworkStats}}.
}
\usage{
GGMnetworkStats.fused(Plist)
}
\arguments{
  \item{Plist}{
  A \code{list} of sparse precision/partial correlation matrix.
}
}
\value{
A \code{data.frame} of the various network statistics for each class. The names of \code{Plist} is prefixed to column-names.
}
\author{
  Anders E. Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\details{
For details on the columns see \code{\link{GGMnetworkStats}}.
}
\seealso{
\code{\link{GGMnetworkStats}}
}
\examples{
## Create some "high-dimensional" data
set.seed(1)
p <- 10
ns <- c(5, 6)
Slist <- createS(ns, p)

## Obtain sparsified partial correlation matrix
Plist    <- ridgeP.fused(Slist, ns, lambda = c(5.2, 1.3), verbose = FALSE)
PCsparse <- sparsify.fused(Plist , threshold = "absValue", absValueCut = 0.2)
SPlist <- lapply(PCsparse, "[[", "sparsePrecision") # Get sparse precisions

## Calculate GGM network statistics in each class
\dontrun{GGMnetworkStats.fused(SPlist)}
}
