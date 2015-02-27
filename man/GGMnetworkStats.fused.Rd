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
Anders Ellern Bilgrau
}
\details{
For details on the columns see \code{\link{GGMnetworkStats}}.
}
\seealso{
\code{\link{GGMnetworkStats}}
}
\examples{
## Create some (high-dimensional) data
p <- 15
ns <- c(10, 11)
Slist <- createS(ns, p)

## Obtain sparsified partial correlation matrix
Plist    <- ridgeS.fused(Slist, ns, lambda = 10, lambdaF = 10, verbose = FALSE)
PCsparse <- sparsify.fused(Plist , threshold = "top")

## Calculate GGM network statistics in each class
GGMnetworkStats.fused(Plist)
}
