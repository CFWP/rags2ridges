\name{GGMblockTest}
\alias{GGMblockTest}
\title{
Test for block-indepedence
}
\description{
Function performing a test that evaluates the null hypothesis of block-independence against the alternative
of block-dependence (presence of non-zero elements in the off-diagonal block) in the precision matrix using high-dimensional data.
The mentioned test is a permutation-based test (see details).
}
\usage{
GGMblockTest(Y, id, nPerm = 1000, lambda, target = default.target(covML(Y)),
             type = "Alt", lowCiThres = 0.1, ncpus = 1, verbose = TRUE)
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
  \item{lambda}{
A \code{numeric} representing the penalty parameter employed in the permutation test.
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{lowCiThres}{
A \code{numeric} taking a value between 0 and 1. Determines speed of efficient p-value calculation.
}
  \item{ncpus}{
A \code{numeric} or \code{integer} indicating the desired number of cpus to be used.
}
  \item{verbose}{
A \code{logical} indicating if information on progress and output should be printed on screen.
}
}
\details{
The function performs a permutation test for the null hypothesis of block-independence against the alternative of
block-dependence (presence of non-zero elements in the off-diagonal block) in the precision matrix using high-dimensional data.
In the low-dimensional setting the common test statistic under multivariate normality (cf. Anderson, 2003) is:
\deqn{
\log( \| \hat{\mathbf{\Sigma}}_a \| ) + \log( \| \hat{\mathbf{\Sigma}}_b \| ) - \log( \| \hat{\mathbf{\Sigma}} \| ),
}
where the \eqn{\hat{\mathbf{\Sigma}}_a}, \eqn{\hat{\mathbf{\Sigma}}_b}, \eqn{\hat{\mathbf{\Sigma}}} are
the estimates of the covariance matrix in the sub- and whole group(s), respectively.

To accommodate the high-dimensionality the parameters of interest are estimated in a penalized manner (ridge-type
penalization, see \code{\link{ridgeP}}). Penalization involves a degree of freedom (the penalty parameter: \code{lambda})
which needs to be fixed before testing. To decide on the penalty parameter for testing we refer to the
\code{\link{GGMblockNullPenalty}} function. With an informed choice on the penalty parameter at hand,
the null hypothesis is evaluated by permutation. Hereto the samples are permutated within block.
The resulting permuted data set represents the null hypothesis. Many permuted data sets are generated.
For each permutation the test statistic is calculated. The observed test statistic is
compared to the null distribution from the permutations.

The function implements an efficient permutation resampling algorithm (see van Wieringen et al., 2008, for details.):
If the probability of a p-value being below \code{lowCiThres} is smaller than 0.001 (read: the test is
unlikely to become significant),
the permutation analysis is terminated and a p-value of unity (1) is reported.

When \code{verbose = TRUE} also graphical output is generated: A histogram of the null-distribution. Note that, when
\code{ncpus} is larger than 1, functionalities from \href{http://cran.r-project.org/package=snowfall}{snowfall}
are imported.
}
\value{
An object of class list:
\item{statistic}{A \code{numeric} representing the observed test statistic (i.e., likelihood ratio).}
\item{pvalue}{A \code{numeric} giving the p-value for the block-independence test.}
\item{nulldist}{A \code{numeric} vector representing the permutation null distribution for the test statistic.}
\item{nperm}{A \code{numeric} indicating the number of permutations used for p-value calculation.}
\item{remark}{A \code{"character"} that states whether the permutation algorithm was terminated prematurely or not.}
}
\references{
Anderson, T.W. (2003). An Introduction to Multivariate Statistical Analysis, 3rd Edition. John Wiley.

van Wieringen, W.N., van de Wiel, M.A., and van der Vaart, A.W. (2008). A Test for Partial Differential Expression.
Journal of the American Statistical Association 103: 1039-1049.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{ridgeP}}, \code{\link{optPenalty.LOOCVauto}}, \code{\link{default.target}}, \code{\link{GGMblockNullPenalty}}
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

## Perform test
testRes <- GGMblockTest(X, id, nPerm = 100, lambdaNull)
}
