% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rags2ridgesFused.R
\name{hist.ptest}
\alias{hist.ptest}
\alias{plot.ptest}
\title{Plot the results of a fusion test}
\usage{
\method{hist}{ptest}(x, add.extra = TRUE, ...)

\method{plot}{ptest}(x, add.extra = TRUE, ...)
}
\arguments{
\item{x}{A \code{ptest} object (a list). Usually the output of
\code{\link{fused.test}}.}

\item{add.extra}{A logical. Add extra information to the plot.}

\item{\dots}{Arguments passed to \code{plot}.}
}
\value{
Invisibly returns \code{x} with extra additions.
}
\description{
Plot a histogram of the null distribution and the observed test statistic in
a permutation type "fusion test".
}
\details{
\code{plot.ptest} is simply a wrapper for \code{hist.ptest}.
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
\references{
Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
  van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
  Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
  of Machine Learning Research, 21(26): 1-52.
}
\seealso{
\code{\link{fused.test}}, \code{\link{print.ptest}}
}
\author{
Anders Ellern Bilgrau, Carel F.W. Peeters <carel.peeters@wur.nl>,
  Wessel N. van Wieringen
}
