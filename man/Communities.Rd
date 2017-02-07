\name{Communities}
\alias{Communities}
\title{
Search and visualize community-structures
}
\description{
Function that searches for and visualizes community-structures in graphs.
}
\usage{
Communities(P, graph = TRUE, lay = "layout_with_fr", coords = NULL,
            Vsize = 15, Vcex = 1, Vcolor = "orangered",
            VBcolor = "darkred", VLcolor = "black", main = "")
}
\arguments{
  \item{P}{
Sparsified precision \code{matrix}
}
  \item{graph}{
A \code{logical} indicating if the results should be visualized.
}
  \item{lay}{
A \code{character} mimicking a call to \code{\link{igraph}} layout functions.
Determines the placement of vertices.
}
  \item{coords}{
A \code{matrix} containing coordinates.
Alternative to the lay-argument for determining the placement of vertices.
}
  \item{Vsize}{
A \code{numeric} determining the vertex size.
}
  \item{Vcex}{
A \code{numeric} determining the size of the vertex labels.
}
  \item{Vcolor}{
A \code{character} (scalar or vector) determining the vertex color.
}
  \item{VBcolor}{
A \code{character} determining the color of the vertex border.
}
  \item{VLcolor}{
A \code{character} determining the color of the vertex labels.
}
  \item{main}{
A \code{character} giving the main figure title.
}
}
\details{
Communities in a network are groups of vertices (modules) that are densely connected within.
Community search is performed by the Girvan-Newman algorithm (Newman and Girvan, 2004).

When \code{graph = TRUE} the community structure in the graph is visualized.
The default layout is according to the Fruchterman-Reingold algorithm (1991).
Most layout functions supported by \code{\link{igraph}} are supported (the function is partly a wrapper around certain \code{\link{igraph}} functions).
The igraph layouts can be invoked by a \code{character} that mimicks a call to a \code{\link{igraph}} layout functions in the \code{lay} argument.
When using \code{lay = NULL} one can specify the placement of vertices with the \code{coords} argument.
The row dimension of this matrix should equal the number of vertices.
The column dimension then should equal 2 (for 2D layouts) or 3 (for 3D layouts).
The \code{coords} argument can also be viewed as a convenience argument as it enables one, e.g., to layout a graph
according to the coordinates of a previous call to \code{Ugraph}.
If both the the lay and the coords arguments are not \code{NULL}, the lay argument takes precedence.
Communities are indicated by color markings.
}
\value{
An object of class list:
\item{membership}{\code{numeric} vector indicating, for each vertex, community membership.}
\item{modularityscore}{\code{numeric} scalar indicating the modularity value of the community structure.}

When \code{graph = TRUE} the function also returns a graph.
}
\references{
Csardi, G. and Nepusz, T. (2006).
The igraph software package for complex network research.
InterJournal, Complex Systems 1695.
http://igraph.sf.net

Fruchterman, T.M.J., and Reingold, E.M. (1991).
Graph Drawing by Force-Directed Placement.
Software: Practice & Experience, 21: 1129-1164.

Newman, M. and Girvan, M. (2004).
Finding and evaluating community structure in networks.
Physical Review E, 69: 026113.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{Ugraph}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)

## Determine support regularized standardized precision under optimal penalty
PC0 <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor

## Search and visualize communities
Commy <- Communities(PC0)
}
