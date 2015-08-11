\name{NLL}
\alias{NLL}
\alias{PNLL}
\alias{NLL.fused}
\alias{PNLL.fused}
\title{
Evaulate the (penalized) (fused) likelihood
}
\description{
Functions that evaulates the (penalized) (fused) likelihood.
}
\usage{
NLL(S, P)
PNLL(S, P, T, lambda)
NLL.fused(Slist, Plist, ns)
PNLL.fused(Slist, Plist, ns, Tlist, lambda)
}
\arguments{
\item{S,Slist}{
A (list of) positive semi definite sample covariance matrices.
}
\item{P,Plist}{
A (list of) positive definite precision matrices.
}
\item{T,Tlist}{
A (list of) positive definite target matrices.
}
\item{ns}{
A \code{numeric} of sample sizes.
}
\item{lambda}{
A \code{numeric} penalty parameter. For the \code{.fused} functions, this is a penalty \code{matrix}.
}
}
\value{
A single number.
}
\author{
Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeP}}, \code{\link{ridgeP.fused}}
}
\examples{
ns <- c(4,5)
Slist <- createS(n = ns, p = 5)
Plist <- list(diag(5), diag(2,5))
Tlist <- list(diag(5), diag(5))

NLL(Slist[[1]], Plist[[1]])
PNLL(Slist[[1]], Plist[[1]], Tlist[[1]], lambda = 1)
NLL.fused(Slist, Plist, ns)
PNLL.fused(Slist, Plist, ns, Tlist, lambda = diag(2))
}
