\name{covML}
\alias{covML}
\title{
Maximum likelihood estimation of the covariance matrix
}
\description{
Function that gives the maximum likelihood estimate of the covariance matrix.
}
\usage{
covML(Y)
}
\arguments{
  \item{Y}{
Data \code{matrix}. Variables assumed to be represented by columns.
}
}
\details{
The function gives the maximum likelihood estimate of the covariance matrix. The input matrix \code{Y} assumes that the variables
are represented by the columns.
}
\value{
Function returns the maximum likelihood estimate of the covariance \code{matrix}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeS}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain ML estimate covariance matrix
Cx <- covML(X)
}