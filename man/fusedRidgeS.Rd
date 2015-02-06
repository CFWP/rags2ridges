\name{fusedRidgeS}
\alias{fusedRidgeS}
\title{
  Fused ridge precision estimation
}
\description{
  Performs the fused ridge estimation of multiple precision matrices in
  cases where multiple classes of data is present for given penalties.
}
\usage{
fusedRidgeS(SList, ns, TList = lapply(SList, default.target),
            lambda1, LambdaP, lambda2,
            PList, maxit = 100L, verbose = TRUE, eps = 1e-04)
}
\arguments{
  \item{SList}{
    A \code{list} of length \eqn{K} of covariance matrices, i.e. square,
    symmetric \code{numeric} matrices of the same size.
    The \eqn{i}th matrix correspond to the \eqn{i}th class.
}
  \item{ns}{
    A \code{numeric} vector of sample sizes on which the matrices in \code{SList}
    are based.
}
  \item{TList}{
    A \code{list} of length \eqn{K} of \code{numeric} target matrices the same sizes
    as in \code{SList}.
    If not supplied, the default is given by \code{\link{default.target}}.
}
  \item{lambda1}{
    A single \code{numeric} giving the regular ridge penalty used.
}
  \item{LambdaP}{
    The fused penalty matrix.
    A symmetric \code{numeric} \code{matrix} of size \eqn{K} times \eqn{K}
    giving the pair-specific penalties.
    I.e. \code{LambdaP[i, j]} is the fused penalty given applied on the pairing
    between classes \eqn{i} and \eqn{j}.
}
  \item{lambda2}{
    A single \code{numeric} giving the fused penalty assuming the complete
    penalty graph, i.e. \code{LambdaP[i,j]} equals \code{lambda2} for all
    \eqn{i} and \eqn{j}.
    Provides an alternative interface to \code{LambdaP}. Either \code{LambdaP}
    or \code{lambda2} must be given.
}
  \item{PList}{
    An optional \code{list} of initial starting values for the algorithm,
    i.e. a list of precision matrices of the same size as \code{SList}.
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
Returns a \code{list} as \code{SList} with precision estimates of the corresponding classes.
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
\code{\link{ridgeS}} for the regular ridge estimate
}
\examples{
# Create some (not at all high-dimensional) data on three classes
p <- 5  # Dimension
ns <- c(4, 6, 8)  # Sample sizes (K = 3 classes)
SList <- createS(ns, p = p)
str(SList, max.level = 2)  # The structure of SList

#
# Estimate the precisions (using the complete penalty graph)
#
res1 <- fusedRidgeS(SList, ns, lambda1 = 0.3, lambda2 = 2)
print(res1)

# The same using the penalty matrix (the diagnal is ignored)
myLambdaP <- matrix(c(0,2,2,
                      2,0,2,
                      2,2,0), 3, 3, byrow = TRUE)
res2 <- fusedRidgeS(SList, ns, lambda1 = 0.3, LambdaP = myLambdaP)
stopifnot(all.equal(res1, res2))



#
# Estimate the precisions (using a non-complete penalty graph)
#

# Say we only want to shrink pairs (1,2) and (2,3) and not (1,3)
myLambdaP[1,3] <- myLambdaP[3,1] <- 0
print(myLambdaP)
res3 <- fusedRidgeS(SList, ns, lambda1 = 0.3, LambdaP = myLambdaP)
# which is not the same as res2 and res1 now!


#
# Using other and custom target matrices
#

# Construct a custom target list
myTList <- list(diag(p), matrix(1, p, p), matrix(0, p, p))
res4 <- fusedRidgeS(SList, ns, TList = myTList, lambda1 = 0.3, lambda2 = 0.2)
print(res4)

# Alternative, see ?default.target
myTList2 <- lapply(SList, default.target, type = "Null")  # For the null target
res5 <- fusedRidgeS(SList, ns, TList = myTList2, lambda1 = 0.3, lambda2 = 0.2)
print(res5)

}

