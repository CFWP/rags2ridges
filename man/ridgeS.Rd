\name{ridgeS}
\alias{ridgeS}
\title{
Ridge estimation for high-dimensional precision matrices
}
\description{
This function is now deprecated. Please use \code{ridgeP} instead.
}
\usage{
ridgeS(S, lambda, type = "Alt", target = default.target(S))
}
\arguments{
  \item{S}{
Sample covariance \code{matrix}.
}
  \item{lambda}{
A \code{numeric} representing the value of the penalty parameter.
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
}
\details{
See \code{ridgeP}.
}
\value{
Function returns a regularized precision \code{matrix}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeP}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Obtain regularized precision matrix
ridgeS(Cx, lambda = 10, type = "Alt")
}
