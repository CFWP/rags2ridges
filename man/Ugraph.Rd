\name{Ugraph}
\alias{Ugraph}
\title{
Visualize undirected graph
}
\description{
Function that visualizes the sparsified precision matrix as an undirected graph.
}
\usage{
Ugraph(M, type = c("plain", "fancy", "weighted"), lay = layout.circle, Vsize = 15,
Vcex = 1, Vcolor = "orangered", VBcolor = "darkred", VLcolor = "black",
prune = FALSE, legend = FALSE, label = "", Lcex = 1.3, PTcex = 4, cut = 0.5,
scale = 10, pEcolor = "black", nEcolor = "grey", main = "")
}
\arguments{
  \item{M}{
(Possibly sparsified) precision \code{matrix}
}
  \item{type}{
A \code{character} indicating the type of graph to be produced. Must be one of: "plain", "fancy", "weighted".
}
  \item{lay}{
Function call to \code{\link{igraph}} determining the placement of vertices.
}
  \item{Vsize}{
A \code{numeric} determining the vertex size.
}
  \item{Vcex}{
A \code{numeric} determining the size of the vertex labels.
}
  \item{Vcolor}{
A \code{character} determining the vertex color.
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
A \code{numeric} representing a scale factor for visualizing strenght of edges when \code{type = "weighted"}.
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
When \code{type = "plain"} a plain undirected graph is given representing the conditional (in)dependencies exemplified by
the sparsified precision.

When \code{type = "fancy"} a more elaborate graph is given in which dashed lines indicate
negative partial correlations while solid lines indicate positive partial correlations, and in which black lines indicate
strong edges. Strong edges are deemed such by setting \code{cut}. If a the absolute value of a precision element \eqn{\geq} \code{cut}
the corresponding edge is deemed strong and colored black in the graph. The argument \code{cut} is thus only used when
\code{type = "fancy"}.

When \code{type = "weighted"} an undirected graph is given in which edge thickness represents the strenght of the partial
correlations. The \code{nEcolor} colored edges then represent negative partial correlations while \code{pEcolor} colored edges
represent positive partial correlations. (Relative) edge thickness in this type of graph can be set by the argument \code{scale}.
The arguments \code{scale}, \code{nEcolor}, and \code{pEcolor} are thus only used when \code{type = "weighted"}.

The default layout gives a circular placement of the vertices. All layout functions supported by
\code{\link{igraph}} are supported (the function is partly a wrapper around certain \code{\link{igraph}} functions).
The legend allows one to specify the kind of variable the vertices represent, such as, e.g., mRNA transcripts. The arguments
\code{label}, \code{Lcex}, and \code{PTcex} are only used when \code{legend = TRUE}.

If \code{prune = TRUE} the vertices of degree 0 (vertices not implicated by any edge) are removed. For the colors supported
by the arguments \code{Vcolor}, \code{VBcolor}, \code{VLcolor}, \code{pEcolor}, and \code{nEcolor} see
\url{www.stat.columbia.edu/~tzheng/files/Rcolor.pdf}.
}
\references{
Csardi, G. and Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems 1695.
http://igraph.sf.net

van Wieringen, W.N. & Peeters, C.F.W. (2015).
Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
arXiv:1403.0904v3 [stat.ME].

van Wieringen, W.N. & Peeters, C.F.W. (2015).
Application of a New Ridge Estimator of the Inverse Covariance Matrix
to the Reconstruction of Gene-Gene Interaction Networks.
In: di Serio, C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.)
`Computational Intelligence Methods for Bioinformatics and Biostatistics'.
Lecture Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
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
}
