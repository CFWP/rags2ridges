\name{ridgeP.fused}
\alias{ridgeP.fused}
\title{
  Fused ridge estimation
}
\description{
  Performs fused ridge estimation of multiple precision matrices in
  cases where multiple classes of data is present for given a penalty matrix.
}
\usage{
ridgeP.fused(Slist, ns, Tlist = default.target.fused(Slist, ns),
             lambda, Plist, maxit = 100L, verbose = TRUE,
             eps = sqrt(.Machine$double.eps))
}
\arguments{
  \item{Slist}{
    A \code{list} of length \eqn{G} of covariance matrices, i.e. square,
    symmetric \code{numeric} matrices of the same size.
    The \eqn{g}th matrix should correspond to the \eqn{g}th class.
  }
  \item{ns}{
    A \code{numeric} vector of sample sizes on which the matrices in
    \code{Slist} are based.
    I.e. \code{ns[g]} correspond to \code{Slist[[g]]}.
  }
  \item{Tlist}{
    A \code{list} of length \eqn{G} of \code{numeric} p.d. target matrices
    corresponding to the matrices in \code{Slist}.
    If not supplied, the default is given by \code{\link{default.target}}.
  }
  \item{lambda}{
    The \eqn{G} by \eqn{G} penalty matrix.
    That is, a symmetric, non-negative \code{numeric} \code{matrix} of size
    \eqn{G} times \eqn{G} giving the class- and pair-specific penalties.
    The diagonal entries are the class specific ridge penalties.
    I.e. \code{lambda[i, i]} is the ridge penalty for class \eqn{i}.
    The off-diagonal entries are the pair-specific fusion penalties.
    I.e. \code{lambda[i, j]} is the fusion penalty applied on the pair
    of classes \eqn{i} and \eqn{j}.

    Alternatively, can be supplied as a \code{numeric} of length 1 or 2.
    If a single number, a diagonal penalty with lambda in the diagonal is used
    If supplied as a \code{numeric} vector of two numbers, the first is used as
    a common ridge penalty and the second as a common fusion penalty.
  }
  \item{Plist}{
    An optional \code{list} of initial precision matrices for the fused ridge
    algorithm the same size as \code{Slist}.
    Can be omitted.
    Default is the nonfused ridge precision estimate using
    the pooled covariance matrix corresponding to setting all fusion penalties
    to zero.
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
    A single positive \code{numeric} giving the convergence threshold.
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
\note{
  For extreme fusion penalties in \code{lambda} the algorithm is quite sensitive to
  the initial values given in \code{Plist}.
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

res1 <- ridgeP.fused(Slist, ns, lambda = c(1.3, 2.1))
print(res1)

# The same using the penalty matrix (the diagnal is ignored)
mylambda  <- matrix(c(1.3, 2.1, 2.1,
                      2.1, 1.3, 2.1,
                      2.1, 2.1, 1.3), 3, 3, byrow = TRUE)
res2 <- ridgeP.fused(Slist, ns, lambda = mylambda)
stopifnot(all.equal(res1, res2))



#
# Estimate the precisions (using a non-complete penalty graph)
#

# Say we only want to shrink pairs (1,2) and (2,3) and not (1,3)
mylambda[1,3] <- mylambda[3,1] <- 0
print(mylambda)
res3 <- ridgeP.fused(Slist, ns, lambda = mylambda)
# which similar to, but not the same as res1 and res2.


#
# Using other custom target matrices
#

# Construct a custom target list
myTlist <- list(diag(p), matrix(1, p, p), matrix(0, p, p))
res4 <- ridgeP.fused(Slist, ns, Tlist = myTlist, lambda = c(1.3, 2.1))
print(res4)

# Alternative, see ?default.target.fused
myTlist2 <- default.target.fused(Slist, ns, type = "Null")  # For the null target
res5 <- ridgeP.fused(Slist, ns, Tlist = myTlist2, lambda = c(1.3, 2.1))
print(res5)
}

