\name{optPenalty.kCVauto}
\alias{optPenalty.kCVauto}
\title{
Automatic search for optimal penalty parameter
}
\description{
Function that performs an 'automatic' search for the optimal penalty parameter for the \code{\link{ridgeP}}
call by employing Brent's method to the calculation of a cross-validated negative log-likelihood score.
}
\usage{
optPenalty.kCVauto(Y, lambdaMin, lambdaMax,
                   lambdaInit = (lambdaMin + lambdaMax)/2,
                   fold = nrow(Y), cor = FALSE,
                   target = default.target(covML(Y)),
                   type = "Alt")
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
  \item{fold}{
A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation.
}
  \item{cor}{
A \code{logical} indicating if the evaluation of the LOOCV score should be performed on the correlation scale.
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
\eqn{K}-fold cross-validated negative log-likelihood score (using a regularized ridge estimator for the precision matrix).
The search for the optimal value is automatic in the sense that in order to invoke the root-finding abilities of the Brent method,
only a minimum value and a maximum value for the penalty parameter need to be specified as well as a starting penalty value.
The value at which the \eqn{K}-fold cross-validated negative log-likelihood score is minimized is deemed optimal.
The function employs the Brent algorithm as implemented in the
\href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html}{optim} function.
}
\value{
An object of class \code{list}:
\item{optLambda}{A \code{numeric} giving the optimal value for the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the precision matrix of the chosen type (see \code{\link{ridgeP}}) under the optimal value of the penalty parameter.}
}
\references{
Brent, R.P. (1971). An Algorithm with Guaranteed Convergence for Finding a Zero of a Function. Computer Journal 14: 422-425.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
When \code{cor = TRUE} correlation matrices are used in the computation of the (cross-validated) negative
log-likelihood score, i.e., the \eqn{K}-fold sample covariance matrix is a matrix on the correlation scale.
When performing evaluation on the correlation scale the data are assumed to be standardized.
If \code{cor = TRUE} and one wishes to used the default target specification one may consider using \code{target = default.target(covML(Y, cor = TRUE))}. This gives a default target under the assumption of standardized data.

Under the default setting of the fold-argument, \code{fold = nrow(Y)}, one performes leave-one-out cross-validation.
}
\seealso{
\code{\link{GGMblockNullPenalty}}, \code{\link{GGMblockTest}}, \code{\link{ridgeP}}, \code{\link{optPenalty.aLOOCV}},
\code{\link{optPenalty.kCV}}, \cr
\code{\link{default.target}}, \code{\link{covML}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty using K = n
OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30); OPT
OPT$optLambda # Optimal penalty
OPT$optPrec   # Regularized precision under optimal penalty

## Another example with standardized data
X <- scale(X, center = TRUE, scale = TRUE)
OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30, cor = TRUE,
                          target = default.target(covML(X, cor = TRUE))); OPT
OPT$optLambda # Optimal penalty
OPT$optPrec   # Regularized precision under optimal penalty

## Another example using K = 5
OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30, fold = 5); OPT
OPT$optLambda # Optimal penalty
OPT$optPrec   # Regularized precision under optimal penalty
}
