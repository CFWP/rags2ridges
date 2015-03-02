\name{ridgeP.fused}
\alias{ridgeP.fused}
\title{
  Fused ridge precision estimation
}
\description{
  Performs the fused ridge estimation of multiple precision matrices in
  cases where multiple classes of data is present for given penalties.
}
\usage{
ridgeP.fused(Slist, ns, Tlist = default.target.fused(Slist, ns),
             lambda, lambdaFmat, lambdaF,
             Plist, maxit = 100L, verbose = TRUE, eps = 1e-04)
}
\arguments{
  \item{Slist}{
    A \code{list} of length \eqn{K} of covariance matrices, i.e. square,
    symmetric \code{numeric} matrices of the same size.
    The \eqn{i}th matrix correspond to the \eqn{i}th class.
}
  \item{ns}{
    A \code{numeric} vector of sample sizes on which the matrices in \code{Slist}
    are based.
}
  \item{Tlist}{
    A \code{list} of length \eqn{K} of \code{numeric} target matrices the same sizes
    as in \code{Slist}.
    If not supplied, the default is given by \code{\link{default.target}}.
}
  \item{lambda}{
    A single \code{numeric} giving the regular ridge penalty used.
}
  \item{lambdaFmat}{
    The fused penalty matrix.
    A symmetric \code{numeric} \code{matrix} of size \eqn{K} times \eqn{K}
    giving the pair-specific penalties.
    I.e. \code{mylambdaFmat[i, j]} is the fused penalty given applied on the pairing
    between classes \eqn{i} and \eqn{j}.
}
  \item{lambdaF}{
    A single \code{numeric} giving the fused penalty assuming the complete
    penalty graph, i.e. \code{lambdaFmat[i,j]} equals \code{lambdaF} for all
    \eqn{i} and \eqn{j}.
    Provides an alternative interface to \code{lambdaFmat}. Either
    \code{lambdaFmat} or \code{lambdaF} must be given.
}
  \item{Plist}{
    An optional \code{list} of initial starting values for the algorithm,
    i.e. a list of precision matrices of the same size as \code{Slist}.
    Can be omitted. Default is the (regular) ridge precision estimate using
    the pooled covariance matrix.
}
  \item{maxit}{
    A single \code{integer} giving the maximum number of allowed iterations.
    Can be set to \code{Inf}.
    If \code{maxit} is hit, a warning is given.
}
  \item{verbose}{
    \code{logical}. Set to \code{TRUE} for extra output.
}
  \item{eps}{
    A single postive \code{numeric} giving the convergence threshold.
}
}
\details{
Performs a coordinate ascent to find the maximum likelihood of the fused
likelihood problem for a given ridge penalty \eqn{lambda} and fused penalty matrix
\eqn{Lambda_f}.
}
\value{
Returns a \code{list} as \code{Slist} with precision estimates of the corresponding classes.
}
\references{
Bilgrau, AE; Peeters CFW; Eriksen, PS; Boegsted, M; & van Wieringen, WN
(in preparation) "Fused Ridge Estimation of Multiple Inverse Covariance
Matrices from High-Dimensional Data Classes"
}
\author{
Anders Ellern Bilgrau \cr
Carel F.W. Peeters <cf.peeters@vumc.nl> \cr
Wessel N. van Wieringen
}
\seealso{
\code{\link{default.penalty}} \cr
\code{\link{ridgeP}} for the regular ridge estimate
}
\examples{
# Create some (not at all high-dimensional) data on three classes
p <- 5  # Dimension
ns <- c(4, 6, 8)  # Sample sizes (K = 3 classes)
Slist <- createS(ns, p = p)
str(Slist, max.level = 2)  # The structure of Slist

#
# Estimate the precisions (using the complete penalty graph)
#
res1 <- ridgeP.fused(Slist, ns, lambda = 0.3, lambdaF = 2)
print(res1)

# The same using the penalty matrix (the diagnal is ignored)
mylambdaFmat <- matrix(c(0,2,2,
                         2,0,2,
                         2,2,0), 3, 3, byrow = TRUE)
res2 <- ridgeP.fused(Slist, ns, lambda = 0.3, lambdaFmat = mylambdaFmat)
stopifnot(all.equal(res1, res2))



#
# Estimate the precisions (using a non-complete penalty graph)
#

# Say we only want to shrink pairs (1,2) and (2,3) and not (1,3)
mylambdaFmat[1,3] <- mylambdaFmat[3,1] <- 0
print(mylambdaFmat)
res3 <- ridgeP.fused(Slist, ns, lambda = 0.3, lambdaFmat = mylambdaFmat)
# which is not the same as res2 and res1 now!


#
# Using other and custom target matrices
#

# Construct a custom target list
myTlist <- list(diag(p), matrix(1, p, p), matrix(0, p, p))
res4 <- ridgeP.fused(Slist, ns, Tlist = myTlist, lambda = 0.3, lambdaF = 0.2)
print(res4)

# Alternative, see ?default.target
myTlist2 <- lapply(Slist, default.target, type = "Null")  # For the null target
res5 <- ridgeP.fused(Slist, ns, Tlist = myTlist2, lambda = 0.3, lambdaF = 0.2)
print(res5)

}

