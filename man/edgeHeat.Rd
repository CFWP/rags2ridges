\name{edgeHeat}
\alias{edgeHeat}
\title{
Visualize (precision) matrix as a heatmap
}
\description{
Function that visualizes a (precision) matrix as a heatmap. May be used to assess visually the elements of a single (possibly sparsified precision)
matrix. May also be used in assessing the performance of edge selection techniques.
}
\usage{
edgeHeat(M, lowColor = "blue", highColor = "red", textsize = 10, diag = TRUE,
         legend = TRUE, main = "")
}
\arguments{
  \item{M}{
(Possibly sparsified precision) \code{matrix}.
}
  \item{lowColor}{
A \code{character} that determines the color scale in the negative range.
}
  \item{highColor}{
A \code{character} that determines the color scale in the positive range.
}
  \item{textsize}{
A \code{numeric} scaling the text size of row and column labels.
}
  \item{diag}{
A \code{logical} determining if the diagonal elements of the matrix should be included in the color scaling. This argument
is only used when \code{M} is a square \code{matrix}.
}
  \item{legend}{
A \code{logical} indicating whether a color legend should be included.
}
  \item{main}{
A \code{character} giving the main figure title.
}
}
\details{
This function utilizes \href{https://cran.r-project.org/package=ggplot2}{ggplot2} (Wickham, 2009) to visualize a matrix as a heatmap: a false color plot in which the
individual matrix entries are represented by colors. \code{lowColor} determines the color scale for matrix entries in the negative range.
\code{highColor} determines the color scale for matrix entries in the positive range. For the colors supported
by the arguments \code{lowColor} and \code{highColor}, see \url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}. White entries in the plot
represent the midscale value of 0. One can opt to set the diagonal entries to the midscale color of white when one is interested in
(heatmapping) the off-diagonal elements only. To achieve this, set \code{diag = FALSE}. Naturally, the \code{diag} argument is only
used when the input matrix \code{M} is a square matrix.

The intended use of the function is to visualize a, possibly sparsified, precision matrix as a heatmap. The function may also be
used, in a graphical modeling setting, to assess the performance of edge selection techniques. However, the function is quite general,
in the sense that it can represent any \code{matrix} as a heatmap.
}
\references{
Wickham, H. (2009). ggplot2: elegant graphics for data analysis. New York: Springer.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{covML}}, \code{\link{ridgeP}}, \code{\link{sparsify}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Obtain regularized precision matrix
P <- ridgeP(Cx, lambda = 10, type = "Alt")

## Obtain sparsified partial correlation matrix
PC0 <- sparsify(P, threshold = "localFDR", FDRcut = .8)$sparseParCor

## Visualize sparsified partial correlation matrix as heatmap
edgeHeat(PC0)
}
