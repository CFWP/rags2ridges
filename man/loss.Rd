\name{loss}
\alias{loss}
\title{
Evaluate regularized precision under various loss functions
}
\description{
Function that evaluates an estimated and possibly regularized precision matrix under various loss functions.
The loss functions are formulated in precision terms. This function may be used to estimate the risk
(vis-a-vis, say, the true precision matrix) of the various ridge estimators employed.
}
\usage{
loss(E, T, precision = TRUE, type = c("frobenius", "quadratic"))
}
\arguments{
  \item{E}{
Estimated (possibly regularized) precision \code{matrix}.
}
  \item{T}{
True (population) covariance or precision \code{matrix}.
}
  \item{precision}{
A \code{logical} indicating if T is a precision matrix.
}
  \item{type}{
A \code{character} indicating which loss function is to be used. Must be one of: "frobenius", "quadratic".
}
}
\details{
Let \eqn{\mathbf{\Omega}} denote a generic \eqn{(p \times p)} population precision matrix and let
\eqn{\hat{\mathbf{\Omega}}(\lambda)} denote a generic ridge estimator of the precision matrix under
generic regularization parameter \eqn{\lambda} (see also \code{\link{ridgeP}}). The function then
considers the following loss functions:
\enumerate{
  \item Squared Frobenius loss, given by:
  \deqn{
    L_{F}[\hat{\mathbf{\Omega}}(\lambda), \mathbf{\Omega}] = \|\hat{\mathbf{\Omega}}(\lambda) -
    \mathbf{\Omega}\|_{F}^{2};
    }
  \item Quadratic loss, given by:
    \deqn{
    L_{Q}[\hat{\mathbf{\Omega}}(\lambda), \mathbf{\Omega}] = \|\hat{\mathbf{\Omega}}(\lambda)
    \mathbf{\Omega}^{-1} - \mathbf{I}_{p}\|_{F}^{2}.
    }
}
The argument \code{T} is considered to be the true precision matrix when \code{precision = TRUE}.
If \code{precision} \code{= FALSE} the argument \code{T} is considered to represent the true covariance matrix.
This statement is needed so that the loss is properly evaluated over the precision, i.e., depending
on the value of the \code{logical} argument \code{precision} inversions are employed where needed.

The function can be employed to assess the risk of a certain ridge precision estimator (see also \code{\link{ridgeP}}).
The risk \eqn{\mathcal{R}_{f}} of the estimator \eqn{\hat{\mathbf{\Omega}}(\lambda)} given a loss function \eqn{L_{f}},
with \eqn{f \in \{F, Q\}} can be defined as the expected loss:
\deqn{
\mathcal{R}_{f}[\hat{\mathbf{\Omega}}(\lambda)] =
\mathrm{E}\{L_{f}[\hat{\mathbf{\Omega}}(\lambda),
\mathbf{\Omega}]\},
}
which can be approximated by the mean or median of losses over repeated simulation runs.
}
\value{
Function returns a \code{numeric} representing the loss under the chosen loss function.
}
\references{
van Wieringen, W.N. & Peeters, C.F.W. (2015).
Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
arXiv:1403.0904v3 [stat.ME].
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{covML}}, \code{\link{ridgeP}}
}
\examples{
## Define population covariance
set.seed(333)
p = 25
n = 1000
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Truecov <- covML(X)

## Obtain sample
samples <- X[sample(nrow(X), 10), ]
Cxx <- covML(samples)

## Obtain regularized precision
P <- ridgeP(Cxx, 10, type = "Alt")

## Evaluate estimated precision against population
## precision under Frobenius loss
loss(P, Truecov, precision = FALSE, type = "frobenius")
}
