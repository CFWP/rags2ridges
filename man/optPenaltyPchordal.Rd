\name{optPenaltyPchordal}
\alias{optPenaltyPchordal}
\title{
Automatic search for penalty parameter of ridge precision estimator with known chordal support
}
\description{
Automic search for the optimal ridge penalty parameter for the ridge estimator of the precision matrix with known chordal support. Optimal in the sense that it yields the maximum cross-validated likelihood. The search employs the Brent algorithm as implemented in the \code{\link[stats:optim]{optim}} function. 
}
\usage{
optPenaltyPchordal(Y, lambdaMin, lambdaMax, 
lambdaInit=(lambdaMin+lambdaMax)/2, 
zeros, cliques=list(), separators=list(), 
target=default.target(covML(Y)), type="Alt")
}
\arguments{
\item{Y}{ Data \code{matrix}. Variables assumed to be represented by columns. }
\item{lambdaMin}{ A \code{numeric} giving the minimum value for the penalty parameter. }
\item{lambdaMax}{ A \code{numeric} giving the maximum value for the penalty parameter. }
\item{lambdaInit}{ A \code{numeric} giving the initial value for the penalty parameter. }
\item{target}{ A target \code{matrix} (in precision terms) for Type I ridge estimators. }
\item{zeros}{ A \code{matrix} with indices of entries of the precision matrix that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the precision matrix. The first column contains the row indices and the second the column indices. The specified conditional independence graph implied by the zero-structure of the precision should be undirected and decomposable. If not, it is symmetrized and triangulated. }
\item{cliques}{ A \code{list}-object containing the node indices per clique as obtained from the \code{\link{support4ridgeP}}-function. }
\item{separators}{ A \code{list}-object containing the node indices per separator as obtained from the \code{\link{support4ridgeP}}-function. }
\item{type}{ A \code{character} indicating the type of ridge estimator to be used. Must be one of: \code{Alt}, \code{ArchI}, \code{ArchII}. }
}
\details{
See the function \code{\link[stats:optim]{optim}} for details on the implementation of the Brent algorithm.
}
\value{
A \code{numeric} with the LOOCV optimal choice for the ridge penalty parameter.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 

Van Wieringen, W.N. and Peeters, C.F.W. (2016), "Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data", \emph{Computational Statistics and Data Analysis}, 103, 284-303.
}
\author{
Wessel N. van Wieringen.
}
\seealso{
\code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}, \code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link[rags2ridges:optPenalty.aLOOCV]{optPenalty.aLOOCV}}, \code{\link[rags2ridges:optPenalty.kCV]{optPenalty.kCV}}
}
\examples{
# generate data
p <- 8
n <- 100
set.seed(333)
Y <- matrix(rnorm(n*p), nrow = n, ncol = p)

# define zero structure  
S <- covML(Y)
S[1:3, 6:8] <- 0
S[6:8, 1:3] <- 0
zeros <- which(S==0, arr.ind=TRUE)

# obtain (triangulated) support info 
supportP <- support4ridgeP(nNodes=p, zeros=zeros)

# determine optimal penalty parameter
\dontrun{
optLambda <- optPenaltyPchordal(Y, 10^(-10), 10, 0.1, zeros=supportP$zeros, 
	cliques=supportP$cliques, separators=supportP$separators)
}
optLambda <- 0.1

# estimate precision matrix with known (triangulated) support
Phat <- ridgePchordal(S, optLambda, zeros=supportP$zeros, 
	cliques=supportP$cliques, separators=supportP$separators)
}
