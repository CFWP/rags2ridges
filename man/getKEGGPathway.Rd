\name{getKEGGPathway}
\alias{getKEGGPathway}
\title{
Download KEGG pathway
}
\description{
Download information and graph object of a given pathway from the Kyoto
Encyclopedia of Genes and Genomes (KEGG) database.
}
\usage{
getKEGGPathway(kegg.id)
}
\arguments{
  \item{kegg.id}{
    A \code{character} giving the KEGG ID, e.g. \code{"map04210"},
    \code{"map04064"}, \code{"map04115"}. Can be prefixed with \code{"path:"}.
  }
}
\value{
  Returns a \code{list} with entries:
  \item{df}{A \code{data.frame} description of the KEGG pathway.}
  \item{graph}{The KEGG pathway represented as a \code{graphNEL} object.}
}
\references{
  \url{https://www.genome.jp/kegg/}
}
\note{
  It is currently necessary to \code{require("KEGGgraph")} (or
  \code{require("KEGGgraph")}) due to a bug in \pkg{KEGGgraph}.
}
\author{
  Anders Ellern Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\details{
  Usage of this function requires an internet connection.
  The igraph objects can be obtained with \code{igraph::igraph.from.graphNEL}.
  The moral graph can be obtained with \code{gRbase::moralize}. To obtain the
  adjacency matrix, use \code{gRbase::as.adjMAT} or
  \code{igraph::get.adjacency}
}
\seealso{
  \code{\link{kegg.target}}
}
\examples{
\dontrun{
if (require("KEGGgraph")) {
  getKEGGPathway("map04064")
}
}
}
\keyword{internal}

