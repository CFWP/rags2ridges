\name{optPenalty.LOOCV}
\alias{optPenalty.LOOCV}
\title{
Select optimal penalty parameter by leave-one-out cross-validation
}
\description{
Function that selects the optimal penalty parameter for the \code{\link{ridgeS}} call by usage of leave-one-out cross-validation.
Its output includes (a.o.) the precision matrix under the optimal value of the penalty parameter.
}
\usage{
optPenalty.LOOCV(Y, lambdaMin, lambdaMax, step, type = "Alt", 
target = default.target(covML(Y)), output = "light", graph = TRUE, 
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
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
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
The function calculates a cross-validated negative log-likelihood score (using a regularized ridge estimator for the precision matrix) for each 
value of the penalty parameter contained in the search grid by way of leave-one-out cross-validation. The value of the penalty parameter 
that achieves the lowest cross-validated negative log-likelihood score is deemed optimal. The penalty parameter must be positive such that \code{lambdaMin} must 
be a positive scalar. The maximum allowable value of \code{lambdaMax} depends on the type of ridge estimator employed. For details on the type of 
ridge estimator one may use (one of: "Alt", "ArchI", "ArchII") see \code{\link{ridgeS}}. 
The ouput consists of an object of class list (see below). When \code{output = "light"} (default) only 
the \code{optLambda} and \code{optPrec} elements of the list are given.
}
\value{
An object of class list:
\item{optLambda}{A \code{numeric} giving the optimal value of the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the precision matrix of the chosen type (see \code{\link{ridgeS}}) under the optimal value of the penalty parameter.}
\item{lambdas}{A \code{numeric} vector representing all values of the penalty parameter for which cross-validation was performed; Only given when \code{output = "all"}.}
\item{LLs}{A \code{numeric} vector representing the mean of cross-validated negative log-likelihoods for each value of the penalty parameter given in \code{lambdas}; Only given when \code{output = "all"}.}
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeS}}, \code{\link{optPenalty.aLOOCV}}, \code{\link{default.target}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT  <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100); OPT	
OPT$optLambda	# Optimal penalty
OPT$optPrec	# Regularized precision under optimal penalty
}