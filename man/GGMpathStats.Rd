\name{GGMpathStats}
\alias{GGMpathStats}
\title{
Gaussian graphical model node pair path statistics
}
\description{
Function that calculates, for a specified node pair representing endpoints, path statistics from a sparse precision matrix.
The sparse precision matrix is taken to represent the conditional independence graph of a Gaussian
graphical model. The contribution to the observed covariance between the specified endpoints is calculated
for each (heuristically) determined path between the endpoints.
}
\usage{
GGMpathStats(P0, node1, node2, neiExpansions = 2, verbose = TRUE, graph = TRUE,
             nrPaths = 2, lay = layout.circle, nodecol = "skyblue", Vsize = 15,
             Vcex = 0.6, VBcolor = "darkblue", VLcolor = "black",
             all.edges = TRUE, prune = TRUE, legend = TRUE, scale = 1,
             Lcex = 0.8, PTcex = 2, main = "")
}
\arguments{
  \item{P0}{
Sparse (possibly standardized) precision matrix.
}
  \item{node1}{
A \code{numeric} specifying an endpoint. The numeric should correspond to a row/column of the precision matrix
and as such represents the corresponding variable.
}
  \item{node2}{
A \code{numeric} specifying a second endpoint. The numeric should correspond to a row/column of the precision matrix
and as such represents the corresponding variable.
}
  \item{neiExpansions}{
A \code{numeric} determining how many times the neighborhood around the respective endpoints should be expanded
in the search for shortest paths between the node pair.
}
  \item{verbose}{
A \code{logical} indicating if a summary of the results should be printed on screen.
}
  \item{graph}{
A \code{logical} indicating if the strongest paths should be visualized with a graph.
}
  \item{nrPaths}{
A \code{numeric} indicating the number of paths (with the highest contribution to the marginal covariance
between the indicated node pair) to be visualized/highlighted.
}
  \item{lay}{
Function call to \code{\link{igraph}} determining the placement of vertices.
}
  \item{nodecol}{
A \code{character} determining the color of \code{node1} and \code{node2}.
}
  \item{Vsize}{
A \code{numeric} determining the vertex size.
}
  \item{Vcex}{
A \code{numeric} determining the size of the vertex labels.
}
  \item{VBcolor}{
A \code{character} determining the color of the vertex borders.
}
  \item{VLcolor}{
A \code{character} determining the color of the vertex labels.
}
  \item{all.edges}{
A \code{logical} indicating if edges other than those implied by the \code{nrPaths}-paths between
\code{node1} and {node2} should also be visualized.
}
  \item{prune}{
A \code{logical} determining if vertices of degree 0 should be removed.
}
  \item{legend}{
A \code{logical} indicating if the graph should come with a legend.
}
  \item{scale}{
A \code{numeric} representing a scale factor for visualizing strenght of edges.
It is a relative scaling factor, in the sense that the edges implied by the \code{nrPaths}-paths between
\code{node1} and {node2} have edge thickness that is twice this scaling factor (so it is
a scaling factor vis-a-vis the unimplied edges).
}
  \item{Lcex}{
A \code{numeric} determining the size of the legend box.
}
  \item{PTcex}{
A \code{numeric} determining the size of the exemplary lines in the legend box.
}
  \item{main}{
A \code{character} giving the main figure title.
}
}
\details{
The conditional independence graph (as implied by the sparse precision matrix) is undirected. In undirected
graphs origin and destination are interchangeable and are both referred to as 'endpoints' of a path. The
function searches for shortest paths between the specified endpoints \code{node1} and \code{node2}.
It searches for shortest paths that visit nodes only once. The shortest paths
between the provided endpoints are determined heuristically by the following procedure. The search is initiated
by application of the \code{get.all.shortest.paths}-function from the \code{\link{igraph}}-package,
which yields all shortest paths between the nodes. Next, the neighborhoods of the endpoints are defined
(excluding the endpoints themselves). Then, the shortest paths are found between: (a)
\code{node1} and node \emph{Vs} in its neighborhood; (b) node \emph{Vs} in the \code{node1}-neighborhood and node
\emph{Ve} in the \code{node2}-neighborhood; and (c) node \emph{Ve} in the \code{node2}-neighborhood and \code{node2}.
These paths are glued and new shortest path candidates are obtained (preserving only novel paths). In additional
iterations (specified by \code{neiExpansions}) the \code{node1}- and \code{node2}-neighborhood are expanded by
including their neighbors (still excluding the endpoints) and shortest paths are again
searched as described above.

The contribution of a particular path to the observed covariance between the specified
node pair is calculated in accordance with Theorem 1 of Jones and West (2005). As in Jones and West (2005),
paths whose weights have an opposite sign to the marginal covariance
(between endnodes of the path) are referred to as 'moderating paths' while paths whose weights
have the same sign as the marginal covariance are referred to as 'mediating' paths. Such
paths are visualized when \code{graph = TRUE}.

All arguments following the \code{graph} argument are only (potentially) used when \code{graph = TRUE}.
When \code{graph = TRUE} the conditional independence graph is returned with the paths highlighted that have the
highest contribution to the marginal covariance between the specified endpoints. The number of paths highlighted
is indicated by \code{nrPaths}. The edges of mediating paths are represented in green while the edges of moderating
paths are represented in red. When \code{all.edges = TRUE} the edges other than those implied by the \code{nrPaths}-paths between
\code{node1} and {node2} are also visualized (in lightgrey). When \code{all.edges = FALSE} only the mediating and
moderating paths implied by \code{nrPaths} are visualized.

The default layout gives a circular placement of the vertices. All layout functions supported by
\code{\link{igraph}} are supported. The arguments \code{Lcex} and \code{PTcex} are only used when \code{legend = TRUE}.
If \code{prune = TRUE} the vertices of degree 0 (vertices not implicated by any edge) are removed. For the colors supported
by the arguments \code{nodecol}, \code{Vcolor}, and \code{VBcolor}, see
\url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}.
}
\value{
An object of class list:
\item{pathStats}{A \code{matrix} specifying the paths, their respective lengths, and their respective contributions to
the marginal covariance between the endpoints.}
\item{paths}{A \code{list} representing the respective paths as numeric vectors.}
\item{Identifier}{A \code{data.frame} in which each numeric from \code{paths} is connected to an identifier such as a variable name.}
}
\note{
Eppstein (1998) describes a more sophisticated algorithm for finding the top \emph{k} shortest paths in a graph.
}
\references{
Eppstein, D. (1998). Finding the k Shortest Paths. SIAM Journal on computing 28: 652-673.

Jones, B., and West, M. (2005). Covariance Decomposition in Undirected Gaussian Graphical Models. Biometrika 92: 779-786.
}
\author{
Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{ridgeP}}, \code{\link{optPenalty.LOOCVauto}}, \code{\link{sparsify}}
}
\examples{
## Obtain some (high-dimensional) data
p <- 25
n <- 10
set.seed(333)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- letters[1:p]

## Obtain regularized precision under optimal penalty
OPT <- optPenalty.LOOCVauto(X, lambdaMin = .5, lambdaMax = 30)

## Determine support regularized standardized precision under optimal penalty
PC0 <- sparsify(OPT$optPrec, threshold = "localFDR")$sparseParCor

## Obtain information on mediating and moderating paths between nodes 14 and 23
pathStats <- GGMpathStats(PC0, 14, 23, verbose = TRUE, prune = FALSE)
pathStats
}
