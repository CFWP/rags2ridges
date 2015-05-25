\name{isSymmetricPD}
\alias{isSymmetricPD}
\title{
  Test for symmetric positive definiteness
}
\description{
  Function to test if a \code{matrix} is symmetric postive definite or not.
}
\usage{
  isSymmetricPD(M)
}
\arguments{
  \item{M}{
   A square symmetric matrix.
  }
}
\value{
  Returns a \code{logical} value. Returns \code{TRUE} if the \code{M} is
  symmetric postive definite and \code{FALSE} if not.
  If \code{M} is not even symmetric, the function throws an error.
}
\details{
  Test positive definiteness by Cholesky decompostion.
}
\author{
  Anders Ellern Bilgrau
  Carel F.W. Peeters <cf.peeters@vumc.nl>, \cr
  Wessel N. van Wieringen
}
\seealso{
  \code{\link{isSymmetric}}
}
\examples{
A <- matrix(rnorm(25), 5, 5)
\dontrun{
isSymmetricPD(A)
}
B <- symm(A)
isSymmetricPD(B)
isSymmetricPD(crossprod(B))
}
