\name{kegg.target}
\alias{kegg.target}
\title{
Construct target matrix from KEGG
}
\description{
Construct a target matrix by combining topology information from the Kyoto Encyclopedia of Genes and Genomes
(KEGG) database and pilot data.
}
\usage{
kegg.target(Y, kegg.id, method = "linreg", organism = "hsa",
            graph = getKEGGPathway(kegg.id)$graph)
}
\arguments{
  \item{Y}{
    The complete observation matrix of observations with variables
    in columns. The column names should be on the form e.g.
    \code{"hsa:3988"} ("\code{<organism>:<Entrez id>}"). It can however also be
    just the Entrez id with or without the post-fixed \code{"_at"} and
    then the specified \code{organism} will be assumed.
  }
  \item{kegg.id}{
    A \code{character} giving the KEGG ID, e.g. \code{"map04210"},
    \code{"map04064"}, or \code{"map04115"}.
  }
  \item{method}{
    The method for estimating the non-zero entries moralized graph of the
    KEGG topology.
    Currently, only \code{"linreg"} is implemented.
  }
  \item{organism}{
    A \code{character} giving the organism, the default is
    \code{"hsa"} (homo-sapiens).
  }
  \item{graph}{
     A \code{graphNEL} object specifying the topology of the pathway.
     Can be used to avoid repeatedly downloading the information.
  }
}
\details{
  The function estimates the precision matrix based on the topology
  given by the KEGG database.
  Requires a connection to the internet.
}
\value{
  Returns a target \code{matrix} with size depending on the \code{kegg.id}.
}
\references{
  \url{http://www.genome.jp/kegg/}
}
\note{
  It is currently nessesary to \code{require("KEGGgraph")} (or
  \code{require("KEGGgraph")}) due to a bug in \pkg{KEGGgraph}.
}
\author{
  Anders Ellern Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\seealso{
  \code{\link{getKEGGPathway}}, \code{\link{default.target}}, and
  \code{\link{default.target.fused}}
}
\examples{
if (require("KEGGgraph")) {
kegg.g <- getKEGGPathway("map04115")$graph

# Create some toy data with the correct names
Y <- createS(n = 10, p = numNodes(kegg.g), dataset = TRUE)
colnames(Y) <- nodes(kegg.g)

T <- kegg.target(Y, "map04115")
print(T[1:10, 1:10])
}
}

