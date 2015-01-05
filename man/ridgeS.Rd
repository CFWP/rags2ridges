\name{ridgeS}
\alias{ridgeS}
\title{
Ridge estimation for high-dimensional precision matrices
}
\description{
Function that calculates various Ridge estimators for high-dimensional precision matrices. 
}
\usage{
ridgeS(S, lambda, type = "Alt", target = default.target(S))
}
\arguments{
  \item{S}{
Sample covariance \code{matrix}.
}
  \item{lambda}{
A \code{numeric} representing the value of the penalty parameter.
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
}
\details{
The function can calculate various ridge estimators for high-dimensional precision matrices. Current (well-known) ridge estimators can be 
roughly divided in two archetypes. The first archetypal form employs a convex combination of \eqn{\mathbf{S}} and a positive definite (p.d.) target matrix
\eqn{\mathbf{T}}: \eqn{\hat{\mathbf{\Omega}}^{\mathrm{I}}(\lambda_{\mathrm{I}}) = [(1-\lambda_{\mathrm{I}}) \mathbf{S} + \lambda_{\mathrm{I}} 
\mathbf{T}]^{-1}}, with \eqn{\lambda_{\mathrm{I}} \in (0,1]}. A common target choice is for \eqn{\mathbf{T}} to be diagonal with
\eqn{(\mathbf{T})_{jj} = (\mathbf{S})_{jj}} for \eqn{j=1, \ldots, p}. The second archetypal form can be given as
\eqn{\hat{\mathbf{\Omega}}^{\mathrm{II}}(\lambda_{\mathrm{II}}) = (\mathbf{S} + \lambda_{\mathrm{II}} \mathbf{I}_{p})^{-1}}
with \eqn{\lambda_{\mathrm{II}} \in (0, \infty)}. Viewed from a penalized estimation perspective, the two archetypes utilize
penalties that do not coincide with the matrix-analogue of the common ridge penalty. van Wieringen and Peeters (2014) derive analytic
expressions for alternative Type I and Type II ridge precision estimators based on a proper L2-penalty. Their alternative Type I estimator
(target shrinkage) takes the form
\deqn{\hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
\left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
\lambda_{a}\mathbf{T})^{2}\right]^{1/2} + \frac{1}{2}(\mathbf{S} -
\lambda_{a}\mathbf{T})\right\}^{-1},}
while their alternative Type II estimator can be given as a special case of the former:
\deqn{\hat{\mathbf{\Omega}}^{\mathrm{II}a}(\lambda_{a}) =
\left\{\left[\lambda_{a}\mathbf{I}_{p} +
\frac{1}{4}\mathbf{S}^{2}\right]^{1/2} + \frac{1}{2}\mathbf{S}\right\}^{-1}.}
These alternative estimators were shown to be superior to the archetypes in terms of risk under various loss functions (van Wieringen and Peeters, 2014).

The \code{lambda} parameter in \code{ridgeS} generically indicates the penalty parameter. It must be chosen in accordance with the type of
ridge estimator employed. The domains for the penalty parameter in the archetypal estimators are given above. The domain for \code{lambda} in
the alternative estimators is \eqn{(0, \infty)}. The \code{type} parameter specifies the type of ridge estimator. Specifying \code{type = "ArchI"} leads
to usage of the archetypal I estimator while specifying \code{type = "ArchII"} leads to usage of the archetypal II estimator. 
In the latter situation the argument \code{target} remains unused. Specifying \code{type = "Alt"}
enables usage of the alternative ridge estimators: when \code{type = "Alt"} and the \code{target} matrix is p.d. one obtains the alternative Type I estimator;
when \code{type = "Alt"} and the \code{target} matrix is specified to be the null-matrix one obtains the alternative Type II estimator.

The Type I estimators thus employ target shrinkage. The default target for both the archetype and alternative is \code{default.target(S)}. When \code{target} is 
not the null-matrix it is expected to be p.d. for the alternative type I estimator. The target is always expected to be p.d. in case of the archetypal I estimator. 
The archetypal Type I ridge estimator is rotation equivariant when the target is of the form \eqn{\mu\mathbf{I}_{p}} with \eqn{\mu \in (0,\infty)}. The archetypal
Type II estimator is rotation equivariant by definition. When the target is of the form \eqn{\varphi\mathbf{I}_{p}} with \eqn{\varphi \in [0,\infty)}, 
then the alternative ridge estimator is rotation equivariant. Its analytic computation is then particularly speedy as the (relatively) expensive matrix square root can then be circumvented. 
}
\value{
Function returns a regularized precision \code{matrix}.
}
\references{
Van Wieringen, W.N. and Peeters, C.F.W. (2014). Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data. arXiv:1403.0904 [stat.ME].
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{default.target}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Obtain regularized precision matrix
ridgeS(Cx, lambda = 10, type = "Alt")
}