\name{rmvnormal}
\alias{rmvnormal}
\title{Multivariate Gaussian simulation}
\description{
Fast simulation from multivariate Gaussian probability distribution.
}
\usage{
rmvnormal(n, mu, sigma)
}
\arguments{
\item{n}{An \code{integer} giving the number of observations to be simulated.}
\item{mu}{A \code{numeric} vector of dimension \eqn{p} giving the means of
  normal distribution.}
\item{sigma}{A variance-covariance \code{matrix} of dimension \eqn{p} times
  \eqn{p}.}
}
\value{
  Returns a \eqn{n} by \eqn{p} matrix of
  observations from a multivariate normal distribution with the given mean
  \code{mu} and covariance
}

\details{
The \code{rmvnormal} function is copied from the \code{GMCM}-package.
It is similar to \code{rmvnorm} from the \code{mvtnorm}-package.
}
\author{
  Anders Ellern Bilgrau
}
\examples{
rmvnormal(n = 10, mu = 1:4, sigma = diag(4))
}

