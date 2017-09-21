\name{KLdiv}
\alias{KLdiv}
\title{
Kullback-Leibler divergence between two multivariate normal distributions
}
\description{
Function calculating the Kullback-Leibler divergence between two multivariate normal distributions.
}
\usage{
KLdiv(Mtest, Mref, Stest, Sref, symmetric = FALSE)
}
\arguments{
  \item{Mtest}{
A \code{numeric} mean vector for the approximating multivariate normal distribution.
}
  \item{Mref}{
A \code{numeric} mean vector for the true/reference multivariate normal distribution.
}
  \item{Stest}{
A covariance \code{matrix} for the approximating multivariate normal distribution.
}
  \item{Sref}{
A covariance \code{matrix} for the true/reference multivariate normal distribution.
}
  \item{symmetric}{
A \code{logical} indicating if the symmetric version of Kullback-Leibler divergence should be calculated.
}
}
\details{
The Kullback-Leibler (KL) information (Kullback and Leibler, 1951; also known as relative entropy) is a measure of divergence between two probability distributions. Typically, one distribution
is taken to represent the `true' distribution and functions as the reference distribution while the other is taken to be an
approximation of the true distribution. The criterion then measures the loss of information in approximating the reference
distribution. The KL divergence between two \eqn{p}-dimensional multivariate normal distributions
\eqn{\mathcal{N}^{0}_{p}(\boldsymbol{\mu}_{0}, \mathbf{\Sigma}_{0})} and \eqn{\mathcal{N}^{1}_{p}(\boldsymbol{\mu}_{1}, \mathbf{\Sigma}_{1})}
is given as
\deqn{
\mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \| \mathcal{N}^{1}_{p}) =
\frac{1}{2}\left\{\mathrm{tr}(\mathbf{\Omega}_{1}\mathbf{\Sigma}_{0})
+ (\boldsymbol{\mu}_{1} - \boldsymbol{\mu}_{0})^{\mathrm{T}}
\mathbf{\Omega}_{1}(\boldsymbol{\mu}_{1} - \boldsymbol{\mu}_{0}) - p
- \ln|\mathbf{\Sigma}_{0}| + \ln|\mathbf{\Sigma}_{1}| \right\},
}
where \eqn{\mathbf{\Omega} = \mathbf{\Sigma}^{-1}}. The KL divergence is not a proper metric as
\eqn{\mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \| \mathcal{N}^{1}_{p}) \neq
\mathrm{I}_{KL}(\mathcal{N}^{1}_{p} \| \mathcal{N}^{0}_{p})}. When \code{symmetric = TRUE} the function calculates
the symmetric KL divergence (also referred to as Jeffreys information), given as
\deqn{
\mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \| \mathcal{N}^{1}_{p}) +
\mathrm{I}_{KL}(\mathcal{N}^{1}_{p} \| \mathcal{N}^{0}_{p}).
}
}
\value{
Function returns a \code{numeric} representing the (symmetric) Kullback-Leibler divergence.
}
\references{
Kullback, S. and Leibler, R.A. (1951). On Information and Sufficiency. Annals of Mathematical Statistics 22: 79-86.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{covML}}, \code{\link{ridgeP}}
}
\examples{
## Define population
set.seed(333)
p = 25
n = 1000
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cov0  <- covML(X)
mean0 <- colMeans(X)

## Obtain sample from population
samples <- X[sample(nrow(X), 10),]
Cov1  <- covML(samples)
mean1 <- colMeans(samples)

## Regularize singular Cov1
P <- ridgeP(Cov1, 10)
CovR <- solve(P)

## Obtain KL divergence
KLdiv(mean1, mean0, CovR, Cov0)
}
