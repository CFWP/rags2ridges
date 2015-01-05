\name{ridgePathS}
\alias{ridgePathS}
\title{
Visualize the regularization path
}
\description{
Function that visualizes the regularization paths of the nonredundant elements of a regularized precision matrix against the 
(range of the) penalty parameter. 
}
\usage{
ridgePathS(S, lambdaMin, lambdaMax, step, type = "Alt", target = default.target(S), 
plotType = "pcor", diag = FALSE, verticle = FALSE, value, verbose = TRUE)
}
\arguments{
  \item{S}{
Sample covariance \code{matrix}.
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
  \item{plotType}{
A \code{character} indicating the type of element for which a visualization of the regularization paths is desired. Must
be one of: "pcor", "cor", "cov", "prec". 
}
  \item{diag}{
A \code{logical} indicating if the diagonal elements should be retained for visualization.
}
  \item{verticle}{
A \code{logical} indicating if output graph should come with a verticle line at a pre-specified value for the penalty parameter.
}
  \item{value}{
A \code{numeric} indicating a pre-specified value for the penalty parameter.
}
  \item{verbose}{
A \code{logical} indicating if information on progress should be printed on screen.
}
}
\details{
The function visualizes the regularization path of the individual elements of a regularized precision matrix against the penalty 
parameter. The range of the penalty parameter is given by [\code{lambdaMin},\code{lambdaMax}]. 
The penalty parameter must be positive such that \code{lambdaMin} must 
be a positive scalar. The maximum allowable value of \code{lambdaMax} depends on the 
type of ridge estimator employed. For details on the type of 
ridge estimator one may use (one of: "Alt", "ArchI", "ArchII") see \code{\link{ridgeS}}. 

Regularization paths may be visualized for (partial) correlations, covariances and precision elements. 
The type of element for which a visualization of the regularization paths is desired can be indicated by the argument \code{plotType}.
When \code{verticle = TRUE} a verticle line is added at the constant \code{value}. This option can be used 
to assess whereabouts the optimal penalty obtained by, e.g., the routines \code{\link{optPenalty.LOOCV}}
or \code{\link{optPenalty.aLOOCV}}, finds itself along the regularization path.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{ridgeS}}, \code{\link{covML}}, \code{\link{optPenalty.LOOCV}}, \code{\link{optPenalty.aLOOCV}}, \code{\link{default.target}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Visualize regularization paths
ridgePathS(Cx, .001, 50, 200, plotType = "pcor")
}