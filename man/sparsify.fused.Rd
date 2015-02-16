\name{sparsify.fused}
\alias{sparsify.fused}

\title{
  Determine support of multiple partial correlation/precision matrices
}
\description{
  A simple wrapper for \code{\link{sparsify}} which determines the support of a \code{list} of partial correlation/precision matrix by various methods and returns the sparsified matrices.
}
\usage{
sparsify.fused(PList, ...)
}
\arguments{
  \item{PList}{
  A \code{list} of \code{numeric} precision matrices.
}
  \item{\dots}{
  Arguments passed to \code{\link{sparsify}}.
}
}
\details{
See \code{\link{sparsify}} for details.
}
\value{
A \code{list} of the same length as \code{PList} with the output from \code{\link{sparsify}}.
}
\author{
Anders Ellern Bilgrau
}
\seealso{
\code{\link{sparsify}}
}
\examples{
ns <- c(10, 11)
YList <- createS(ns, p = 16, covariance = FALSE)
SList <- lapply(YList, covML)
TList <- default.target.fused(SList, ns)

# Obtain regularized precision under optimal penalty
opt <- optPenalty.fused.LOOCVauto(YList, TList, approximate = TRUE,
                                  maxit.ridgeS.fused = 1500)
# Use the optimal penalties
PList <- ridgeS.fused(SList, ns, lambda = opt$lambda,
                      lambdaFmat = opt$lambdaFmat, maxit = 1000)

# Determine support regularized (standardized) precision under optimal penalty
res <- sparsify.fused(PList, threshold = "top", verbose = FALSE)
round(res[[1]]$sparsePrecision, 1)
round(res[[2]]$sparsePrecision, 1)
}

