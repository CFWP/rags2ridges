\name{ridgePsign}
\alias{ridgePsign}
\title{
Ridge estimation for high-dimensional precision matrices with known sign of off-diagonal precision elements.
}
\description{
Function that calculates the ridge estimators for high-dimensional precision matrices with known sign of the off-diagonal precision elements.
}
\usage{
ridgePsign(S, lambda, sign, target=default.target(S), method="nlm", verbose=TRUE, ...)
}
\arguments{
\item{S}{ Sample covariance \code{matrix}. }
\item{lambda}{ A \code{numeric} representing the value of the penalty parameter. }
\item{target}{ A target \code{matrix} (in precision terms) for the ridge precision estimator. }
\item{sign}{ A character indicating the required sign of the off-diagonal elements of ridge precision estimate. Must be either: "pos" (positive) and "neg" (negative). }
\item{method}{ A \code{character} : which optimization function should be used: \code{"nlm"} (default) or \code{"optim"} which refer to \code{\link[stats:nlminb]{nlminb}} or \code{\link[stats:constrOptim]{constrOptim}}, respectively. }
\item{verbose}{ \code{Logical} indicator: should intermediate output be printed on the screen? }
\item{...}{ Additional arguments passed on to either \code{\link[stats:nlminb]{nlminb}} or \code{\link[stats:constrOptim]{constrOptim}}.}
}
\details{
Modified version of the \code{\link{ridgePchordal}}-function, now the ridge precision matrix estimate has off-diagonal elements equalling zero or of the specified sign. The estimate is found by solving a constrained estimation problem. This is done numerically and employs the \code{\link[stats:nlminb]{nlminb}} and \code{\link[stats:constrOptim]{constrOptim}} procedure of R. These procedures are initiated by the ridge ML precision estimate and its off-diagonal elements with the excluded sign set to (effectively) zero.
}
\value{
The function returns a regularized precision \code{matrix} with off-diagonal elements of specified signed or zero.
}
\author{
W.N. van Wieringen.
}
\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}
}
\examples{
# obtain some data
p <- 8
n <- 100
set.seed(333)
Y <- matrix(rnorm(n*p), nrow = n, ncol = p)

# obtain regularized precision matrix with off-diagonal elements of specified signed 
ridgePsign(covML(Y), lambda=0.1, sign="pos")
}
