\name{Ugraph}
\alias{Ugraph}
\title{
Visualize undirected graph
}
\description{
Function that visualizes the sparsified precision matrix as an undirected graph.
}
\usage{
Ugraph(M, type = c("plain", "fancy", "weighted"),
       lay = "layout_in_circle", coords = NULL, Vsize = 15,
       Vcex = 1, Vcolor = "orangered", VBcolor = "darkred",
       VLcolor = "black", prune = FALSE, legend = FALSE,
       label = "", Lcex = 1.3, PTcex = 4, cut = .5,
       scale = 10, pEcolor = "black", nEcolor = "grey",
       main = "")
}
\arguments{
  \item{M}{
(Possibly sparsified) precision \code{matrix}
}
  \item{type}{
A \code{character} indicating the type of graph to be produced. Must be one of: "plain", "fancy", "weighted".
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
  \item{prune}{
A \code{logical} determining if vertices of degree 0 should be removed.
}
  \item{legend}{
A \code{logical} indicating if the graph should come with a legend.
}
    \item{label}{
A \code{character} giving a name to the legend label.
}
    \item{Lcex}{
A \code{numeric} determining the size of the legend box.
}
    \item{PTcex}{
A \code{numeric} determining the size of the exemplary vertex in the legend box.
}
    \item{cut}{
A \code{numeric} indicating the cut-off for indicating strong edges when \code{type = "fancy"}.
}
    \item{scale}{
A \code{numeric} representing a scale factor for visualizing strength of edges when \code{type = "weighted"}.
}
  \item{pEcolor}{
A \code{character} determining the color of the edges tied to positive precision elements. Only when \code{type = "weighted"}.
}
  \item{nEcolor}{
A \code{character} determining the color of the edges tied to negative precision elements. Only when \code{type = "weighted"}.
}
  \item{main}{
A \code{character} giving the main figure title.
}
}
\details{
The intended use of this function is to visualize a sparsified precision/partial correlation matrix as an undirected graph.
When \code{type = "plain"} a plain undirected graph is given representing the conditional (in)dependencies exemplified by the sparsified precision.

When \code{type = "fancy"} a more elaborate graph is given in which dashed lines indicate
negative partial correlations while solid lines indicate positive partial correlations, and in which grey lines indicate strong edges.
Strong edges are deemed such by setting \code{cut}.
If a the absolute value of a precision element \eqn{\geq} \code{cut} the corresponding edge is deemed strong and colored grey in the graph.
The argument \code{cut} is thus only used when \code{type = "fancy"}.

When \code{type = "weighted"} an undirected graph is given in which edge thickness represents the strength of the partial correlations.
The \code{nEcolor} colored edges then represent negative partial correlations while \code{pEcolor} colored edges
represent positive partial correlations.
(Relative) edge thickness in this type of graph can be set by the argument \code{scale}.
The arguments \code{scale}, \code{nEcolor}, and \code{pEcolor} are thus only used when \code{type = "weighted"}.

The default layout gives a circular placement of the vertices. Most layout functions supported by
\code{\link{igraph}} are supported (the function is partly a wrapper around certain \code{\link{igraph}} functions).
The igraph layouts can be invoked by a \code{character} that mimicks a call to a \code{\link{igraph}} layout functions in the \code{lay} argument.
When using \code{lay = NULL} one can specify the placement of vertices with the \code{coords} argument.
The row dimension of this matrix should equal the number of (pruned) vertices.
The column dimension then should equal 2 (for 2D layouts) or 3 (for 3D layouts).
The \code{coords} argument can also be viewed as a convenience argument as it enables one, e.g., to layout a graph
according to the coordinates of a previous call to \code{Ugraph}.
If both the the lay and the coords arguments are not \code{NULL}, the lay argument takes precedence

The legend allows one to specify the kind of variable the vertices represent, such as, e.g., mRNA transcripts.
The arguments \code{label}, \code{Lcex}, and \code{PTcex} are only used when \code{legend = TRUE}.

If \code{prune = TRUE} the vertices of degree 0 (vertices not implicated by any edge) are removed.
For the colors supported by the arguments \code{Vcolor}, \code{VBcolor}, \code{VLcolor}, \code{pEcolor}, and \code{nEcolor} see \url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}.
}
\value{
The function returns a graph.
The function also returns a \code{matrix} object containing the coordinates of the vertices in the given graph.
}
\references{
Csardi, G. and Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems 1695.
http://igraph.sf.net

van Wieringen, W.N. & Peeters, C.F.W. (2016).
Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
Computational Statistics & Data Analysis, vol. 103: 284-303.
Also available as arXiv:1403.0904v3 [stat.ME].

van Wieringen, W.N. & Peeters, C.F.W. (2015).
Application of a New Ridge Estimator of the Inverse Covariance Matrix
to the Reconstruction of Gene-Gene Interaction Networks.
In: di Serio, C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.)
`Computational Intelligence Methods for Bioinformatics and Biostatistics'.
Lecture Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{ridgeP}}, \code{\link{optPenalty.LOOCV}}, \code{\link{optPenalty.aLOOCV}}, \code{\link{sparsify}}
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

## Obtain graphical representation
Ugraph(PC0, type = "fancy", cut = 0.07)

## Obtain graphical representation with Fruchterman-Reingold layout
Ugraph(PC0, type = "fancy", lay = "layout_with_fr", cut = 0.07)

## Add pruning
Ugraph(PC0, type = "fancy", lay = "layout_with_fr",
       cut = 0.07, prune = TRUE)

## Obtain graph and its coordinates
Coordinates <- Ugraph(PC0, type = "fancy", lay = "layout_with_fr",
                      cut = 0.07, prune = TRUE)
Coordinates
}
