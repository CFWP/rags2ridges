\name{sparsify.fused}
\alias{sparsify.fused}

\title{
  Determine support of multiple partial correlation/precision matrices
}
\description{
  A simple wrapper for \code{\link{sparsify}} which determines the support of a \code{list} of partial correlation/precision matrix by various methods and returns the sparsified matrices.
}
\usage{
sparsify.fused(Plist, ...)
}
\arguments{
  \item{Plist}{
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
A \code{list} of the same length as \code{Plist} with the output from \code{\link{sparsify}}.
}
\author{
Anders Ellern Bilgrau,
Wessel N. van Wierigen,
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{sparsify}}
}
\examples{
ns <- c(10, 11)
Ylist <- createS(ns, p = 16, dataset = TRUE)
Slist <- lapply(Ylist, covML)
Tlist <- default.target.fused(Slist, ns)

# Obtain regularized precision under optimal penalty
opt <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "aLOOCV",
                            maxit.ridgeP.fused = 1500)
# Use the optimal penalties
Plist <- ridgeP.fused(Slist, ns, lambda = opt$lambda, maxit = 1000)

# Determine support regularized (standardized) precision under optimal penalty
res <- sparsify.fused(Plist, threshold = "top", verbose = FALSE)
round(res[[1]]$sparsePrecision, 1)
round(res[[2]]$sparsePrecision, 1)
}

