\name{print.ptest}
\alias{print.ptest}
\alias{summary.ptest}
\title{
Print and summarize fusion test
}
\description{
Print and summary functions for the fusion test performed by \code{\link{fused.test}}.
}
\usage{
\method{print}{ptest}(x, digits = 4L, ...)

\method{summary}{ptest}(object, ...)
}
\arguments{
  \item{x, object}{
    The object to print or summarize. Usually the output of \code{\link{fused.test}}.
  }
  \item{digits}{
    An \code{integer} controlling the number of printed digits.
  }
  \item{\dots}{
    Arguments passed on.
    In \code{summary.ptest} the arguments are passed to \code{print.ptest}.
    In \code{print.ptest} are passed to the standard \code{summary} function.
  }
}
\value{
Invisibly returns the object.
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
\code{\link{fused.test}}, \code{\link{hist.ptest}}
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
}
