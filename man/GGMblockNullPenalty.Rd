\name{GGMblockNullPenalty}
\alias{GGMblockNullPenalty}
\title{
Generate the distribution of the penalty parameter under the null hypothesis of block-independence
}
\description{
Function that serves as a precursor function to the block-independence test (see \code{\link{GGMblockTest}}).
It generates an empirical distribution of the penalty parameter under the null hypothesis of block independence 
(in the regularized precision matrix). 
}
\usage{
GGMblockNullPenalty(Y, id, nPerm = 25, lambdaMin, lambdaMax, 
lambdaInit = (lambdaMin + lambdaMax)/2, target = default.target(covML(Y)), 
type = "Alt", ncpus = 1)
}
\arguments{
  \item{Y}{
Data \code{matrix}. Variables assumed to be represented by columns.
}
  \item{id}{
A \code{numeric} vector acting as an indicator variable for two blocks of the precision matrix. The blocks should be coded as \code{0} and \code{1}.
}
  \item{nPerm}{
A \code{numeric} or \code{integer} determining the number of permutations.
}
  \item{lambdaMin}{
A \code{numeric} giving the minimum value for the penalty parameter.
}
  \item{lambdaMax}{
A \code{numeric} giving the maximum value for the penalty parameter.
}
  \item{lambdaInit}{
A \code{numeric} giving the initial value for the penalty parameter for starting optimization.
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{ncpus}{
A \code{numeric} or \code{integer} indicating the desired number of cpus to be used.
}
}
\details{
This function can be viewed as a precursor to the function for the block-independence test (see \code{\link{GGMblockTest}}). 
The mentioned test evaluates the null hypothesis of block-independence against the alternative 
of block-dependence (presence of non-zero elements in the off-diagonal block) in the precision matrix using high-dimensional data. 
To accommodate the high-dimensionality the parameters of interest are estimated in a penalized manner (ridge-type 
penalization, see \code{\link{ridgeS}}). Penalization involves a degree of freedom (the penalty parameter) which needs to 
be fixed before testing. This function then generates an empirical distribution of this penalty parameter. 
Hereto the samples are permutated within block. The resulting permuted data sets represent the null hypothesis. 
To avoid the dependence on a single permutation, many permuted data sets are generated. For each permutation 
the optimal penalty parameter is determined by means of cross-validation (see \code{\link{optPenalty.LOOCVauto}}). 
The resulting optimal penalty parameters are returned. An estimate of the location (such as the median) is recommended for use in the 
block-independence test.
}
\value{
A \code{numeric} vector, representing the distribution of the (LOOCV optimal) penalty parameter under the null 
hypothesis of block-independence. 
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{ridgeS}}, \code{\link{optPenalty.LOOCVauto}}, \code{\link{default.target}}, \code{\link{GGMblockTest}}
}
\examples{
## Obtain some (high-dimensional) data
p = 15
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:15] = letters[1:15]
id <- c(rep(0, 10), rep(1, 5))

## Generate null distribution of the penalty parameter
lambda0dist <- GGMblockNullPenalty(X, id, 5, 0.001, 10)

## Location of null distribution
lambdaNull <- median(lambda0dist)
}
