\name{pooledS}
\alias{pooledS}
\title{
Compute pooled covariance
}
\description{
Compute the pooled covariance estimate from a \code{list} of covariance matrices.
}
\usage{
pooledS(Slist, ns, mle = TRUE, subset = rep(TRUE, length(ns)))
}
\arguments{
  \item{Slist}{A \code{list} of length \eqn{K} of \code{numeric} covariance
    matrices of the same size.}
  \item{ns}{A \code{numeric} vector for length \eqn{K} giving the sample sizes
    in the corresponding entries of \code{Slist}}
  \item{mle}{\code{logical}. If \code{TRUE}, the (biased) MLE is given. If
    \code{FALSE}, the biased corrected estimate is given. Default is
    \code{TRUE}.}
  \item{subset}{\code{logical} vector of the same length as \code{Slist} giving
    the classes to pool. Default is all classes.}
}
\details{
When \code{mle} is \code{FALSE} the given covariance matrices is assumed to have been computed using the denominator \code{ns[i] - 1}. Hence, the sum of all \code{ns} minus \eqn{K} is used a the denominator of the pooled estimate.
Conversely, when \code{mle} is \code{TRUE} the total sum of the sample sizes \code{ns} is used as the denominator in the pooled estimate.
}
\value{
Returns a \code{numeric} matrix with the same size as the elements of \code{Slist}.
}
\author{
Anders Ellern Bilgrau
}
\examples{
ns <- c(4, 6, 8)
Slist <- createS(ns, p = 6)

pooledS(Slist, ns)
pooledS(Slist, ns, mle = FALSE)

# Pool the first two classes only, leave out the remaning
pooledS(Slist, ns, subset = c(TRUE, TRUE, FALSE))
pooledS(Slist, ns, subset = ns > 5) # Pool studies with sample size > 5
}

