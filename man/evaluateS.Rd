\name{evaluateS}
\alias{evaluateS}
\title{
Evaluate numerical properties square matrix
}
\description{
Function that evaluates various numerical properties of a square input matrix. The intended use is to evaluate the various numerical properties of 
what is assumed to be a covariance matrix. Another use is to evaluate the various numerical properties of a (regularized) precision matrix.
}
\usage{
evaluateS(S, verbose = TRUE)
}
\arguments{
  \item{S}{
Covariance or (regularized) precision \code{matrix}.
}
  \item{verbose}{
A \code{logical} indicating if output should be printed on screen.
}
}
\details{
The function evaluates various numerical properties of a covariance or precision input matrix. The function assesses if the input matrix is symmetric,
if all its eigenvalues are real, if all its eigenvalues are strictly positive, and if it is a diagonally dominant matrix. In addition,
the function calculates the trace, the determinant, and the spectral condition number of the input matrix. See, e.g., Harville (1997) 
for more details on the mentioned (numerical) matrix properties.
}
\value{
\item{symm}{A \code{logical} indicating if the matrix is symmetric.}
\item{realEigen}{A \code{logical} indicating if the eigenvalues are real.}
\item{posEigen}{A \code{logical} indicating if the eigenvalues are strictly positive.}
\item{dd}{A \code{logical} indicating if the matrix is diagonally dominant.}
\item{trace}{A \code{numerical} giving the value of the trace.}
\item{det}{A \code{numerical} giving the value of the determinant.}
\item{condNumber}{A \code{numerical} giving the value of the spectral condition number.}
}
\references{
Harville, D.A.(1997). Matrix algebra from a statistician's perspective. New York: Springer-Verlag.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{covML}}, \code{\link{ridgeS}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Evaluate numerical properties covariance matrix
## Obtain, e.g., value trace
Seval <- evaluateS(Cx); Seval
Seval$trace		

## Evaluate numerical properties precision matrix after regularization
P <- ridgeS(Cx, lambda = 10, type = 'Alt')
Peval <- evaluateS(P); Peval
}