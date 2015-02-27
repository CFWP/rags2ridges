\name{default.target.fused}
\alias{default.target.fused}
\title{
  Generate data-driven targets for fused ridge estimation
}
\description{
  Generates a list of (data-driven) targets to use in fused ridge estimation
  Simply a wrapper for \code{\link{default.target}}.
}
\usage{
default.target.fused(Slist, ns, type = "DAIE", equal = TRUE, ...)
}
\arguments{
  \item{Slist}{
    A \code{list} of length \eqn{K} of \code{numeric} covariance matrices of the same size for \eqn{K} classes.
}
  \item{ns}{
    A \code{numeric} vector of sample sizes corresponding to the entries of \code{Slist}.
}
  \item{type}{
    A \code{character} giving the choice of target to construct. See \code{\link{default.target}} for the available options. Default is \code{"DAIE"}.
}
  \item{equal}{
    \code{logical.} If \code{TRUE} the chosen targets are based on the pooled covariance estimate and all targets will be the same. If \code{FALSE}, the \eqn{i}th target is based on the \eqn{i}th entry of \code{Slist}. Default is \code{TRUE}.
}
  \item{\dots}{
    Arguments passed to \code{\link{default.target}}.
}
}
\value{
  A \code{list} of \eqn{K} covariance target matrices of the same size.
}
\author{
  Anders E. Bilgrau, \cr
  Carel F.W. Peeters <cf.peeters@vumc.nl>, \cr
  Wessel N. van Wieringen
}
\seealso{
\code{\link{default.target}}
}
\examples{
# Make some toy data
ns <- c(3, 4)  # Two classes with sample size 3 and 4
Slist <- createS(ns, p = 3)  # Generate two 3-dimensional covariance matrices
Slist

# Different choices:
default.target.fused(Slist, ns)
default.target.fused(Slist, ns, equal = FALSE)
default.target.fused(Slist, ns, type = "Null")
default.target.fused(Slist, ns, type = "DAPV")
default.target.fused(Slist, ns, type = "DAPV", equal = FALSE)
}
