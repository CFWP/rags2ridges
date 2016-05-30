\name{is.Xlist}
\alias{is.Xlist}
\title{
Test if fused list-formats are correctly used
}
\description{
Function to check if the argument submits to the various \code{list}-formats used by the fused ridge estimator and related functions are correct.
That is, it tests if generic fused list arguments (such as \code{Slist}, \code{Tlist}, \code{Plist}, \code{Ylist}) are properly formatted.
}
\usage{
is.Xlist(Xlist, Ylist = FALSE, semi = FALSE)
}
\arguments{
  \item{Xlist}{
    A \code{list} of precision matrices of equal size (\code{Plist}), sample covariance matrices (\code{Slist}), data matrices (\code{Ylist})
  }
  \item{Ylist}{
    \code{logical}. Is \code{Xlist} a \code{list} of data matrices with the same number of columns (\code{Ylist}).
  }
  \item{semi}{
    \code{logical}. Should the matrices in the list be tested to be positive semi definite or positive definite?
  }
}
\value{
  Returns \code{TRUE} if all tests are passed, throws error if not.
}
\references{
  Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and van Wieringen, W.N. (2015).
  Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes,
  arXiv:1509.07982v1 [stat.ME].

  van Wieringen, W.N. & Peeters, C.F.W. (To appear).
  Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
  Computational Statistics & Data Analysis.
  Available as arXiv:1403.0904v3 [stat.ME].
}
\author{
Anders Bilgrau,
Carel F.W. Peeters <cf.peeters@vumc.nl>,
Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeP.fused}}, \code{\link{optPenalty.fused}}
}
\examples{
Slist <- createS(n = c(4, 6, 9), p = 10)
is.Xlist(Slist, semi = TRUE)
}
