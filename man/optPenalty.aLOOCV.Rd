\name{optPenalty.aLOOCV}
\alias{optPenalty.aLOOCV}
\title{
Select optimal penalty parameter by approximate leave-one-out cross-validation
}
\description{
Function that selects the optimal penalty parameter for the \code{\link{ridgeP}} call by usage of approximate
leave-one-out cross-validation. Its output includes (a.o.) the precision matrix under the optimal value of the
penalty parameter.
}
\usage{
optPenalty.aLOOCV(Y, lambdaMin, lambdaMax, step, type = "Alt",
                  cor = FALSE, target = default.target(covML(Y)),
                  output = "light", graph = TRUE, verbose = TRUE)
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
  \item{cor}{
A \code{logical} indicating if the evaluation of the approximate LOOCV score should be performed on the correlation
scale.
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
The function calculates an approximate leave-one-out cross-validated (aLOOCV) negative log-likelihood score (using a
regularized ridge estimator for the precision matrix) for each value of the penalty parameter contained in the search
grid. The utilized aLOOCV score was proposed by Lian (2011) and Vujacic et al. (2014). The aLOOCV negative log-likeliho
od score is computationally more efficient than its non-approximate counterpart (see \code{\link{optPenalty.LOOCV}}).
For details on the aLOOCV negative log-likelihood score see Lian (2011) and Vujacic et al (2014).
For scalar matrix  targets (see \code{\link{default.target}}) the complete solution path of the alternative
Type I and II ridge estimators (see \code{\link{ridgeP}}) depends on only 1 eigendecomposition and 1 matrix
inversion, making the determination of the optimal penalty value particularly efficient (see van Wieringen
and Peeters, 2015).

The value of the penalty parameter that achieves the lowest aLOOCV negative log-likelihood score is
deemed optimal. The penalty parameter must be positive such that \code{lambdaMin} must be a
positive scalar. The maximum allowable value of \code{lambdaMax} depends on the type of
ridge estimator employed. For details on the type of ridge estimator one may use
(one of: "Alt", "ArchI", "ArchII") see \code{\link{ridgeP}}. The ouput consists
of an object of class list (see below). When \code{output = "light"} (default)
only the \code{optLambda} and \code{optPrec} elements of the list are given.
}
\value{
An object of class list:
\item{optLambda}{A \code{numeric} giving the optimal value of the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the precision matrix of the chosen type (see \code{\link{ridgeS}}) under the optimal value of the penalty parameter.}
\item{lambdas}{A \code{numeric} vector representing all values of the penalty parameter for which approximate cross-validation was performed; Only given when \code{output = "all"}.}
\item{aLOOCVs}{A \code{numeric} vector representing the approximate cross-validated negative log-likelihoods for each value of the penalty parameter given in \code{lambdas}; Only given when \code{output = "all"}.}
}
\references{
Lian, H. (2011). Shrinkage tuning parameter selection in precision matrices estimation. Journal of Statistical
Planning and Inference, 141: 2839-2848.

van Wieringen, W.N. & Peeters, C.F.W. (2016).
Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
Computational Statistics & Data Analysis, vol. 103: 284-303.
Also available as arXiv:1403.0904v3 [stat.ME].

Vujacic, I., Abbruzzo, A., and Wit, E.C. (2014). A computationally fast alternative to cross-validation in
penalized Gaussian graphical models. arXiv: 1309.6216v2 [stat.ME].
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\note{
When \code{cor = TRUE} correlation matrices are used in the computation of the approximate (cross-validated) negative
log-likelihood score, i.e., the sample covariance matrix is a matrix on the correlation scale.
When performing evaluation on the correlation scale the data are assumed to be standardized.
If \code{cor = TRUE} and one wishes to used the default target specification one may consider using \code{target = default.target(covML(Y, cor = TRUE))}. This gives a default target under the assumption of standardized data.
}
\seealso{
\code{\link{ridgeP}}, \code{\link{optPenalty.LOOCV}}, \code{\link{optPenalty.LOOCVauto}}, \cr
\code{\link{default.target}}, \code{\link{covML}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT  <- optPenalty.aLOOCV(X, lambdaMin = .001, lambdaMax = 30, step = 400); OPT
OPT$optLambda	# Optimal penalty
OPT$optPrec	  # Regularized precision under optimal penalty

## Another example with standardized data
X <- scale(X, center = TRUE, scale = TRUE)
OPT  <- optPenalty.aLOOCV(X, lambdaMin = .001, lambdaMax = 30,
                          step = 400, cor = TRUE,
                          target = default.target(covML(X, cor = TRUE))); OPT
OPT$optLambda	# Optimal penalty
OPT$optPrec	  # Regularized precision under optimal penalty
}
