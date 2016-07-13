\name{fullMontyS}
\alias{fullMontyS}
\title{
Wrapper function
}
\description{
Function that forms a wrapper around certain \code{rags2ridges} functionalities.
More specifically, it (automatically) invokes functionalities to get from high-dimensional
data to a penalized precision estimate, to the corresponding conditional independence graph
and topology summaries.
}
\usage{
fullMontyS(Y, lambdaMin, lambdaMax, target = default.target(covML(Y)),
           dir = getwd(), fileTypeFig = "pdf", FDRcut = 0.9,
           nOutput = TRUE, verbose = TRUE)
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
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{dir}{
A \code{character} specifying the directory in which the (visual) output is to be stored.
}
  \item{fileTypeFig}{
A \code{character} determining the file type of visual output. Must be one of: "pdf", "eps".
}
  \item{FDRcut}{
A \code{numeric} indicating the cut-off for partial correlation element selection based on local FDR thresholding.
}
  \item{nOutput}{
A \code{logical} indicating if numeric output should be returned.
}
  \item{verbose}{
A \code{logical} indicating if progress updates should be printed on screen.
}
}
\details{
The wrapper always uses the alternative ridge precision estimator (see \code{\link{ridgeP}}) with \code{target}
as the target matrix. The optimal value for the penalty parameter is determined by employing Brent's method to
the calculation of a cross-validated negative log-likelihood score (see \code{\link{optPenalty.LOOCVauto}}).
The support of the regularized precision matrix is determined by way of local FDR thresholding
(see \code{\link{sparsify}}). The corresponding conditional independence graph is visualized using
\code{\link{Ugraph}} with \code{type = "fancy"}. This visualization as well as the calculation
of network statistics (see \code{\link{GGMnetworkStats}}) is based on the standardization of the regularized
and sparsified precision matrix to a partial correlation matrix.
}
\value{
The function stores in the specified directory \code{dir} a condition number plot (either .pdf or .eps file),
a visualization of the network (either .pdf or .eps file), and a file containing network statistics (.txt file).
When \code{nOutput = TRUE} the function also returns an object of class \code{list}:
\item{optLambda}{A \code{numeric} giving the optimal value of the penalty parameter.}
\item{optPrec}{A \code{matrix} representing the regularized precision matrix under the optimal value of the penalty parameter.}
\item{sparseParCor}{A \code{matrix} representing the sparsified partial correlation matrix.}
\item{networkStats}{A \code{matrix} giving the calculated network statistics.}
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\note{
We consider this to be a preliminary version of an envisioned wrapper than will take
better form with subsequent versions of \code{rags2ridges}.
}
\seealso{
\code{\link{ridgeP}}, \code{\link{conditionNumberPlot}}, \code{\link{optPenalty.LOOCVauto}}, \code{\link{sparsify}},
\code{\link{Ugraph}}, \code{\link{GGMnetworkStats}}
}
\examples{
\dontrun{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Employ the wrapper function
theWorks <- fullMontyS(X, lambdaMin = .5, lambdaMax = 30)}
}
