\name{fused.test}
\alias{fused.test}
\title{
Test the necessity of fusion
}
\description{
Function for testing the null hypothesis that all population precision
matrices are equal and thus the necessity for the fusion penalty.
Note, the test performed is conditional on the supplied penalties and targets.
}
\usage{
fused.test(Ylist, Tlist, lambda, n.permutations = 100,
           verbose = FALSE, ...)
}
\arguments{
  \item{Ylist}{
    A \code{list} of length \eqn{G} of observations matrices for each class.
    Variables are assumed to correspond to the columns.
  }
  \item{Tlist}{
    A \code{list} of target matrices for each class. Should be same length
    as \code{Ylist}-
  }
  \item{lambda}{
    A non-negative, symmetric \eqn{G} by \eqn{G} \code{matrix} giving the
    ridge and fusion penalties.
  }
  \item{n.permutations}{
    The number of permutations to approximate the null distribution.
    Default is 100. Should be increased if sufficient computing power is
    available.
  }
  \item{verbose}{
    Print out extra progress information
  }
  \item{\dots}{
    Arguments passed to \code{\link{ridgeP.fused}}.
  }
}
\details{
The function computes the observed score statistic \eqn{U_obs} using the fused
ridge estimator on the given data.
Next, the score statistic is computed a number of times (given by
\code{n.permutations}) under the null-hypothesis by effectively permuting the
class labels of the data.
}
\value{
Returns a \code{list} values containing the observed test statistic and
the test statistic under the null distribution.
}
\references{
  Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and van Wieringen, W.N. (2020).
  Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes.
  Journal of Machine Learning Research, 21(26): 1-52.
}
\author{
Anders Ellern Bilgrau,
Carel F.W. Peeters <cf.peeters@vumc.nl>,
Wessel, N. van Wieringen
}
\seealso{
\code{\link{ridgeP.fused}}
}
\examples{
ns <- c(10, 5, 23)
Ylist <- createS(ns, p = 15, topology = "banded", dataset = TRUE)

# Use the identity target matrix for each class
Tlist <- replicate(length(ns), diag(15), simplify = FALSE)

# Do the test
lm <- matrix(10, 3, 3)
diag(lm) <- 1
ft <- fused.test(Ylist, Tlist, lambda = lm,
                 n.permutations = 500)
print(ft)

# Summary spits out a bit more information
summary(ft)

# The returned object can alo be plotted via
hist(ft)
# or via the alias
plot(ft)

# Customization and parameters work a usual:
hist(ft, col = "steelblue", main = "Null distribution", add.extra = FALSE,
     xlab = "Score statistic", freq = FALSE)
}
