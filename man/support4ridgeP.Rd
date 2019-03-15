\name{support4ridgeP}
\alias{support4ridgeP}
\title{
Support of the adjacency matrix to cliques and separators.
}
\description{
Convert the support of an undirected, chordal graph into a lists of cliques and separators. When the graph is not chordal, it is triangulated to make it so. The undirected graph may be specified as an adjacency matrix, or by the complement of its support as a matrix with the indices of the adjancency matrix corresponding to absent edges. The function thus caters for the two different types of output from the \code{\link[rags2ridges:sparsify]{sparsify}}-function. The function is meant to preceede the \code{\link{ridgePchordal}}, as it its output directly feeds into the latter.
}
\usage{
support4ridgeP(adjMat=NULL, nNodes=NULL, zeros=NULL, verbose=FALSE)
}
\arguments{
\item{adjMat}{ Adjacency matrix of an undirected graph. }
\item{nNodes}{ Positive \code{integer} of length one: number nodes of the network. }
\item{zeros}{ A \code{matrix} with indices of entries of the adjacency matrix that are zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
}
\details{
Essentially, it is a wrapper for the \code{rip}-function from the \code{gRbase}-package, which takes different input and yields slightly different output. Its main purpose is to mold the input such that it is convenient for the \code{ridgePchordal}-function, which provides ridge maximum likelihood estimation of the precision matrix with known support.
}
\value{
A \code{list}-object comprising three slots: 'zeros', 'cliques, 'separators' and 'addedEdges'. The 'zeros'-slot: a \code{matrix} with indices of entries of the adjacency matrix that are zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, it is symmetrized and triangulated. Hence, it may differ from the input 'zeros'. The 'cliques'-slot: a \code{list}-object containing the node indices per clique as obtained from the \code{rip}-function. The 'separators'-slot: a \code{list}-object containing the node indices per clique as obtained from the \code{rip}-function. The 'addedEdges'-slot: a \code{matrix} with indices of edges that have been added in the triangulation.
}
\references{
Lauritzen, S.L. (2004). \emph{Graphical Models}. Oxford University Press.
}
\author{
Wessel N. van Wieringen.
}
\seealso{
\code{\link[rags2ridges:sparsify]{sparsify}},
\code{\link{ridgePchordal}},
\code{\link[gRbase:graph-rip]{gRbase::rip}}.
}
\examples{
# obtain some (high-dimensional) data
p <- 8
n <- 100
set.seed(333)
Y <- matrix(rnorm(n*p), nrow = n, ncol = p)

# create sparse precision
P <- covML(Y)
P[1:3, 6:8] <- 0
P[6:8, 1:3] <- 0

# draw some data
S <- covML(matrix(rnorm(n*p), nrow = n, ncol = p))

# obtain (triangulated) support info
zeros <- which(P==0, arr.ind=TRUE)
supportP <- support4ridgeP(adjMat=adjacentMat(P))

# alternative specification of the support
zeros <- which(P==0, arr.ind=TRUE)
supportP <- support4ridgeP(nNodes=p, zeros=zeros)

# estimate precision matrix with known (triangulated) support
Phat <- ridgePchordal(S, 0.1, zeros=supportP$zeros,
	cliques=supportP$cliques, separators=supportP$separators)
}
