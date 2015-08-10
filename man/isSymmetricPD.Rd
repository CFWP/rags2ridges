\name{isSymmetricPD}
\alias{isSymmetricPD}
\alias{isSymmetricPSD}
\title{
  Test for symmetric positive definiteness
}
\description{
  Function to test if a \code{matrix} is symmetric positive definite or not.
}
\usage{
isSymmetricPD(M)
isSymmetricPSD(M, tol = 1e-4)
}
\arguments{
  \item{M}{
  A square symmetric matrix.
  }
  \item{tol}{
  A numeric giving the tolerance for determining positive semi definiteness.
  }
}
\value{
  Returns a \code{logical} value. Returns \code{TRUE} if the \code{M} is
  symmetric positive (semi)definite and \code{FALSE} if not.
  If \code{M} is not even symmetric, the function throws an error.
}
\details{
  Tests positive definiteness by Cholesky decomposition.
  Tests positive semi definiteness by checking if all eigenvalues are larger than \eqn{-\epsilon|\lambda_1|} where \eqn{\epsilon} is the tolerance and \eqn{\lambda_1} is the largest eigenvalue.
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

C <- crossprod(B)
isSymmetricPD(C)

isSymmetricPSD(C)
}
