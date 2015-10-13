\name{default.target.fused}
\alias{default.target.fused}
\title{
  Generate data-driven targets for fused ridge estimation
}
\description{
  Generates a list of (data-driven) targets to use in fused ridge estimation.
  Simply a wrapper for \code{\link{default.target}}.
}
\usage{
default.target.fused(Slist, ns, type = "DAIE", by, ...)
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
  \item{by}{
    A \code{character} vector with the same length as
    \code{Slist} specifying which groups should share target.
    For each unique entry of \code{by} a target is constructed.
    If omitted, the default is to assign a unique target to each class.
    If not given as a \code{character} coercion into one is attempted.
  }
  \item{\dots}{
    Arguments passed to \code{\link{default.target}}.
}
}
\value{
  A \code{list} of \eqn{K} covariance target matrices of the same size.
}
\author{
  Anders E. Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
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
default.target.fused(Slist, ns, by = seq_along(Slist)) # The same as before
default.target.fused(Slist, ns, type = "Null")
default.target.fused(Slist, ns, type = "DAPV")
default.target.fused(Slist, ns, type = "DAPV", by = rep(1, length(Slist)))


# Make some (more) toy data
ns <- c(3, 4, 6, 7)  # Two classes with sample size 3 and 4
Slist <- createS(ns, p = 2)  # Generate four 2-dimensional covariance matrices

# Use the same target in class 1 and 2, but another in class 3 and 4:
default.target.fused(Slist, ns, by = c("A", "A", "B", "B"))

}
