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
createS(n, p, topology = c("identity", "star", "clique", "chain", "banded"), 
        covariance = TRUE)
}
\arguments{
  \item{n}{A \code{numeric} vector giving number of samples. If the length is larger than 1, the covariance matrices are returned as a list.}
  \item{p}{A \code{numeric} of length 1 giving the dimension of the samples/covariance.}
  \item{topology}{character. The topology to use for the simulations. See the details.}
  \item{covariance}{A \code{logical} value specifying whether the sample covariance or the simulated data itself should be returned.}
}
\details{
The data is simulated from a zero-mean \code{p}-dimensional multivariate gaussian distribution with some precision matrix determined by the argument \code{topology} which defines the GGM.
For the available values of \code{topology} are described below.
\itemize{
\item \code{"identity"} uses the identity matrix (\code{diag(p)}).
\item Others topologies are not implemented yet.
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
# Generate some sample covariance matrices
createS(n = 10, p = 3)
createS(n = c(3, 4, 5), p = 3)
createS(n = c(32, 55), p = 7)

# Generate some (possibly high-dimensional) datasets
createS(n = c(3, 4), p = 6, covariance = FALSE)
createS(n = 3, p = 2, cov = FALSE)
}

