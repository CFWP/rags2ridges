\name{ridgePchordal}
\alias{ridgePchordal}
\title{
Ridge estimation for high-dimensional precision matrices with known chordal support
}
\description{
Function that calculates various ridge estimators for high-dimensional precision matrices with known support. This support should form a chordal graph. If the provided support is not chordal, the function makes it so.
}
\usage{
ridgePchordal(S, lambda, zeros, cliques=list(), separators=list(), 
target=default.target(S), type="Alt", optimizer="nlm", grad=FALSE, 
verbose=TRUE, ...)
}
\arguments{
\item{S}{ Sample covariance \code{matrix}. }
\item{lambda}{ A \code{numeric} representing the value of the penalty parameter. }
\item{target}{ A target \code{matrix} (in precision terms) for Type I ridge estimators. }
\item{zeros}{ \code{Matrix} with indices of entries of the adjacency matrix that are zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, use the \code{\link{support4ridgeP}} to symmetrize and triangulate. This is done automatically if \code{cliques} and \code{separators} arguments are empty lists (and the then employed \code{zeros}-object may differ from the one provided as input). }
\item{cliques}{ A \code{list}-object containing the node indices per clique as obtained from the \code{\link{support4ridgeP}}-function. }
\item{separators}{ A \code{list}-object containing the node indices per separator as obtained from the \code{\link{support4ridgeP}}-function. }
\item{type}{ A \code{character} indicating the type of ridge estimator to be used. Must be one of: \code{Alt} (default), \code{ArchI}, \code{ArchII}. }
\item{optimizer}{ A \code{character} (either \code{nlm} (default) or \code{optim}) specifying which optimization function should be used: \code{\link[stats:nlm]{nlm}} (default) or \code{\link[stats:optim]{optim}}? }
\item{grad}{ A \code{logical} indicator: should, next to the precision matrix estimate, also the gradient be returned? }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{ Additional arguments passed on to either \code{\link[stats:nlm]{nlm}} or \code{\link[stats:optim]{optim}}.}
}
\details{
Sister function to the \code{\link{ridgeP}}-function, incorporating a chordal zero structure of the precision matrix.

The loss function for \code{type="ArchII"} is:
\deqn{ \log(| \mathbf{\Omega} |) - \mbox{tr} ( \mathbf{S} \mathbf{\Omega} ) + \lambda \big\{ \log(| \mathbf{\Omega} |) - \mbox{tr} [ (\mathbf{S} + (1+\lambda) \mathbf{I}_{p \times p}) \mathbf{\Omega} ] \big\}. }
For \code{type="ArchI"} it is:
\deqn{ (1-\lambda) \big [ \log(| \mathbf{\Omega} |) - \mbox{tr} ( \mathbf{S} \mathbf{\Omega} ) \big] + \lambda \big[ \log(| \mathbf{\Omega} |) - \mbox{tr} ( \mathbf{\Omega} ) \big], }
which is obtained from:
\deqn{ \log(| \mathbf{\Omega} |) - \mbox{tr} ( \mathbf{S} \mathbf{\Omega} ) + \nu \big[ \log(| \mathbf{\Omega} |) - \mbox{tr} ( \mathbf{\Omega} ) \big] }
by division of \eqn{(1+\nu)} and writing \eqn{\lambda = \nu / (1 + \nu)}.

An explicit expression for the minimizer of the loss functions implied by the archetypal ridge estimators (\code{type="ArchI"} and \code{type="ArchII"}) exists. For the simple case in which the graph decomposes into cliques 
\eqn{\mathcal{C}_1}, \eqn{\mathcal{C}_2} and separator \eqn{\mathcal{S}} the estimator is:
\deqn{\widehat{\mathbf{\Omega}}  = 
\left(
\begin{array}{lll}
\, [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{C}_1 \setminus \mathcal{S}, \mathcal{C}_1 \setminus \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{C}_1 \setminus \mathcal{S}, \mathcal{S}} & \mathbf{0}_{|\mathcal{C}_1 \setminus \mathcal{S}| \times |\mathcal{C}_2 \setminus \mathcal{S}|}
\\
\, [\widehat{\mathbf{\Omega}}^{(\mathcal{C}_1)}]_{\mathcal{S}, \mathcal{C}_1 \setminus \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{S}, \mathcal{S}} + [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{S}, \mathcal{S}} - \widehat{\mathbf{\Omega}}^{(\mathcal{S})} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{S}, \mathcal{C}_2 \setminus \mathcal{S}}
\\
\, \mathbf{0}_{|\mathcal{C}_2 \setminus \mathcal{S}| \times |\mathcal{C}_1 \setminus \mathcal{S}|} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{C}_2 \setminus \mathcal{S}, \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{C}_2 \setminus \mathcal{S}, \mathcal{C}_2 \setminus \mathcal{S}}
\end{array}
\right),
}
where \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}}, \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}} and \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{S}})}} are the marginal ridge ML covariance estimators for cliques \eqn{\mathcal{C}_1}, \eqn{\mathcal{C}_2} and separator \eqn{\mathcal{S}}. The general form of the estimator, implemented here, is analogous to that provided in Proposition 5.9 of Lauritzen (2004). The proof that this estimator indeed optimizes the corresponding loss function is fully analogous to that of Proposition 5.6 of Lauritzen (2004).  
 
In case, \code{type="Alt"} no explicit expression of the maximizer of the ridge penalized log-likelihood exists. However, an initial estimator analogous to that for \code{type="ArchI"} and \code{type="ArchII"} can be defined. In various boundary cases (\eqn{\lambda=0}, \eqn{\lambda=\infty}, and \eqn{\mathcal{S} = \emptyset}) this initial estimator actually optimizes the loss function. In general, however, it does not. Nevertheless, it functions as well-educated guess for any Newton-like optimization method: convergence is usually achieved quickly. The Newton-like procedure optimizes an unconstrained problem equivalent to that of the penalized log-likelihood with known zeros for the precision matrix (see Dahl \emph{et al}., 2005 for details). 
}

\value{
If \code{grad=FALSE}, the function returns a regularized precision \code{matrix} with specified chordal sparsity structure.

If \code{grad=TRUE}, a list is returned comprising of \emph{i)} the estimated precision matrix, and \emph{ii)} the gradients at the initial and at the optimal (if reached) value. The gradient is returned and it can be checked whether it is indeed (close to) zero at the optimum.
}
\references{
Dahl, J., Roychowdhury, V., Vandenberghe, L. (2005), "Maximum likelihood estimation of Gaussian graphical models: numerical implementation and topology selection", Technical report, UCLA, 2005.

Lauritzen, S.L. (2004). \emph{Graphical Models}. Oxford University Press.  

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191.
}
\author{
Wessel N. van Wieringen.
}
\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}}
}
\examples{
# obtain some (high-dimensional) data
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

# estimate precision matrix with known (triangulated) support
Phat <- ridgePchordal(S, 0.1, zeros=supportP$zeros, 
	cliques=supportP$cliques, separators=supportP$separators)
}
