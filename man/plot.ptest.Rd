\name{plot.ptest}
\alias{plot.ptest}
\alias{hist.ptest}
\title{
  Plot the results of a fusion test
}
\description{
Plot a histogram of the null distribution and the observed test statistic in a permutation type "fusion test".
}
\usage{
## S3 method for class 'ptest'
\method{plot}{ptest}(x, add.extra = TRUE, ...)

## S3 method for class 'ptest'
\method{hist}{ptest}(x, add.extra = TRUE, ...)
}
\arguments{
  \item{x}{
    A \code{ptest} object (a list). Usually the output of \code{\link{fused.test}}.
  }
  \item{add.extra}{
    A logical. Add extra information to the plot.
  }
  \item{\dots}{
    Arguments passed to \code{plot}.
  }
}
\details{
  \code{plot.ptest} is simply a wrapper for \code{hist.ptest}.
}
\value{
  Invisibly returns \code{x} with extra additions.
}
\references{
  Bilgrau, AE; Peeters CFW; Eriksen, PS; Boegsted, M; & van Wieringen, WN
  (In preparation) "Targeted Fused Ridge Estimation of Multiple Inverse Covariance
  Matrices from High-Dimensional Data Classes"
}
\author{
  Anders Ellern Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\seealso{
\code{\link{fused.test}}, \code{\link{print.ptest}}
}
\examples{
ns <- c(10, 5, 23)
Ylist <- createS(ns, p = 15, topology = "banded", dataset = TRUE)

# Use the identity target matrix for each class
Tlist <- replicate(length(ns), diag(15), simplify = FALSE)

# Do the test
lam <- matrix(10, 3, 3)
diag(lam) <- 1
ft <- fused.test(Ylist, Tlist, lambda = lam, n.permutations = 500)

# The returned object can alo be plotted via
hist(ft)
# or via the alias
plot(ft)
}

