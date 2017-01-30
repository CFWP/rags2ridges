\name{optPenalty.kCV}
\alias{optPenalty.kCV}
\title{
Select optimal penalty parameter by \eqn{K}-fold cross-validation
}
\description{
Function that selects the optimal penalty parameter for the \code{\link{ridgeP}} call
by usage of \eqn{K}-fold cross-validation.
Its output includes (a.o.) the precision matrix under the optimal value of the penalty parameter.
}
\usage{
optPenalty.kCV(Y, lambdaMin, lambdaMax, step, fold = nrow(Y),
               cor = FALSE, target = default.target(covML(Y)),
               type = "Alt", output = "light", graph = TRUE,
               verbose = TRUE)
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
  \item{step}{
An \code{integer} determining the number of steps in moving through the grid [\code{lambdaMin}, \code{lambdaMax}].
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
  \item{output}{
A \code{character} indicating if the output is either heavy or light. Must be one of: "all", "light".
}
  \item{graph}{
A \code{logical} indicating if the grid search for the optimal penalty parameter should be visualized.
}
  \item{verbose}{
A \code{logical} indicating if information on progress should be printed on screen.
}
}
\details{
The function calculates a cross-validated negative log-likelihood score (using a regularized ridge
estimator for the precision matrix) for each
value of the penalty parameter contained in the search grid by way of
\eqn{K}-fold cross-validation. The value of the penalty parameter
that achieves the lowest cross-validated negative log-likelihood score
is deemed optimal. The penalty parameter must be positive such that \code{lambdaMin} must
be a positive scalar. The maximum allowable value of \code{lambdaMax}
depends on the type of ridge estimator employed. For details on the type of
ridge estimator one may use (one of: "Alt", "ArchI", "ArchII") see \code{\link{ridgeP}}.
The ouput consists of an object of class list (see below). When \code{output = "light"} (default) only
the \code{optLambda} and \code{optPrec} elements of the list are given.
}
\value{
An object of class list:
\item{optLambda}{A \code{numeric} giving the optimal value of the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the precision matrix of the chosen type (see \code{\link{ridgeP}}) under the optimal value of the penalty parameter.}
\item{lambdas}{A \code{numeric} vector representing all values of the penalty parameter for which cross-validation was performed; Only given when \code{output = "all"}.}
\item{LLs}{A \code{numeric} vector representing the mean of cross-validated negative log-likelihoods for each value of the penalty parameter given in \code{lambdas}; Only given when \code{output = "all"}.}
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\note{
When \code{cor = TRUE} correlation matrices are used in the computation of the (cross-validated) negative
log-likelihood score, i.e., the \eqn{K}-fold sample covariance matrix is a matrix on the correlation scale.
When performing evaluation on the correlation scale the data are assumed to be standardized.
If \code{cor = TRUE} and one wishes to used the default target specification one may consider using \code{target = default.target(covML(Y, cor = TRUE))}. This gives a default target under the assumption of standardized data.

Under the default setting of the fold-argument, \code{fold = nrow(Y)}, one performes leave-one-out cross-validation.
}
\seealso{
\code{\link{ridgeP}}, \code{\link{optPenalty.kCVauto}}, \code{\link{optPenalty.aLOOCV}}, \cr
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
OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100); OPT
OPT$optLambda	# Optimal penalty
OPT$optPrec	  # Regularized precision under optimal penalty

## Another example with standardized data
X <- scale(X, center = TRUE, scale = TRUE)
OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100, cor = TRUE,
                       target = default.target(covML(X, cor = TRUE))); OPT
OPT$optLambda	# Optimal penalty
OPT$optPrec	  # Regularized precision under optimal penalty

## Another example using K = 5
OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100, fold = 5); OPT
OPT$optLambda	# Optimal penalty
OPT$optPrec	  # Regularized precision under optimal penalty
}
