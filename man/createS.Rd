\name{createS}
\alias{createS}
\title{
  Simulate sample covariances or datasets
}
\description{
Simulate data from a p-dimensional (zero-mean) gaussian graphical model (GGM) and return the sample covariance matrix or matrices.
Can also return the original simulated data.
}
\usage{
createS(n, p,
        topology = c("identity", "complete", "star", "clique", "chain",
                     "banded", "barabassi", "small-world"),
        dataset = FALSE, nonzero = 0.25, m = 1L)
}
\arguments{
  \item{n}{A \code{numeric} vector giving number of samples. If the length is larger than 1, the covariance matrices are returned as a list.}
  \item{p}{A \code{numeric} of length 1 giving the dimension of the samples/covariance.}
  \item{topology}{character. The topology to use for the simulations. See the details.}
  \item{dataset}{A \code{logical} value specifying whether the sample covariance or the simulated data itself should be returned.}
  \item{nonzero}{A \code{numeric} of length 1 giving the value of the nonzero entries used in some typologies.}
  \item{m}{A \code{integer} giving the number of blocks (i.e. conditionally independent components) to create. If \code{m} is greater than 1, then the given \code{topology} is used on \code{m} blocks of approximately equal size.}
}
\details{
The data is simulated from a zero-mean \code{p}-dimensional multivariate gaussian distribution with some precision matrix determined by the argument \code{topology} which defines the GGM.
For the available values of \code{topology} are described below.
\itemize{
\item \code{"identity"} uses the identity matrix (\code{diag(p)}) as precision matrix.
  Corresponds to no conditional dependencies.
\item \code{"complete"} is a complete graph with
\item The remaning topologies are not implemented yet.
}
}
\value{
The returned type is dependent on \code{n} and \code{covariance}.
The function generally returns a \code{list} of \code{numeric} matrices with the same length as \code{n}.
If \code{covariance} is \code{FALSE} the simulated datasets with size \code{n[i]} by \code{p} are given in the \code{i} entry of the output.
If \code{covariance} is \code{TRUE} the \code{p} by \code{p} sample covariances of the datasets are given.

When \code{n} has length 1 the \code{list} structure is dropped and the matrix is returned.
}
\author{
Anders Ellern Bilgrau
}
\examples{
## Generate some sample covariance matrices
createS(n = 10, p = 3)
createS(n = c(3, 4, 5), p = 3)
createS(n = c(32, 55), p = 7)

## Generate some datasets and not sample covariance matrices
createS(c(3, 4), p = 6, dataset = TRUE)

## Generate sample covariance matrices from other topologies:
A <- createS(3000, p = 4, topology = "star")
round(solve(A), 3)
B <- createS(3000, p = 4, topology = "banded")
round(solve(B), 3)
C <- createS(3000, p = 4, topology = "clique")  # The complete graph (as m = 1)
round(solve(C), 3)
D <- createS(3000, p = 4, topology = "chain")
round(solve(D), 3)

## Generate smaple covariance matrices from block topologies:
C3 <- createS(3000, p = 10, topology = "clique", m = 3)
round(solve(C3), 1)
C5 <- createS(3000, p = 10, topology = "clique", m = 5)
round(solve(C5), 1)
}

