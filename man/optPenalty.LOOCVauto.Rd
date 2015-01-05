\name{optPenalty.LOOCVauto}
\alias{optPenalty.LOOCVauto}
\title{
Automatic search for optimal penalty parameter
}
\description{
Function that performs an 'automatic' search for the optimal penalty parameter for the \code{\link{ridgeS}}
call by employing Brent's method to the calculation of a cross-validated negative log-likelihood score.
}
\usage{
optPenalty.LOOCVauto(Y, lambdaMin, lambdaMax, lambdaInit = (lambdaMin + lambdaMax)/2, 
target = default.target(covML(Y)), type = "Alt")
}
\arguments{
  \item{Y}{
Data \code{matrix}. Variables assumed to be represented by columns.
}
  \item{lambdaMin}{
A \code{numeric} giving the minimum value for the penalty parameter.
}
  \item{lambdaMax}{
A \code{numeric} giving the maximum value for the penalty parameter.
}
  \item{lambdaInit}{
A \code{numeric} giving the initial (starting) value for the penalty parameter.
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
}
\details{
The function determines the optimal value of the penalty parameter by application of the Brent algorithm (1971) to the 
(leave-one-out) cross-validated negative log-likelihood score (using a regularized ridge estimator for the precision matrix).
The search for the optimal value is automatic in the sense that in order to invoke the root-finding abilities of the Brent method,
only a minimum value and a maximum value for the penalty parameter need to be specified as well as a starting penalty value. 
The value at which the (leave-one-out) cross-validated negative log-likelihood score is minimized is deemed optimal.
The function employs the Brent algorithm as implemented in the 
\href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html}{optim} function. 
}
\value{
An object of class \code{list}:
\item{optLambda}{A \code{numeric} giving the optimal value for the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the precision matrix of the chosen type (see \code{\link{ridgeS}}) under the 
optimal value of the penalty parameter.}
}
\references{
Brent, R.P. (1971). An Algorithm with Guaranteed Convergence for Finding a Zero of a Function. Computer Journal 14: 422-425.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{GGMblockNullPenalty}}, \code{\link{GGMblockTest}}, \code{\link{ridgeS}}, \code{\link{optPenalty.aLOOCV}}, 
\code{\link{optPenalty.LOOCV}}, \code{\link{default.target}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT  <- optPenalty.LOOCVauto(X, lambdaMin = .001, lambdaMax = 30); OPT	
OPT$optLambda # Optimal penalty
OPT$optPrec # Regularized precision under optimal penalty
}
