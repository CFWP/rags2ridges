\name{print.optPenaltyFusedGrid}
\alias{print.optPenaltyFusedGrid}
\alias{plot.optPenaltyFusedGrid}
\title{
Print and plot functions for fused grid-based cross-validation
}
\description{
Print and plot functions for the output from \code{\link{optPenalty.fused.grid}} which performs a grid based cross-validation (CV) search to find optimal penalty parameters.
Currently, only the complete penalty graph is supported.
}
\usage{
\method{print}{optPenaltyFusedGrid}(x, ...)

\method{plot}{optPenaltyFusedGrid}(x, add.text = TRUE, add.contour = TRUE, col = rainbow(100, end = 0.8), ...)
}
\arguments{
  \item{x}{
    A \code{optPenaltyFusedGrid}-object print or plot.
    Usually the output of
    \code{\link{optPenalty.fused.grid}}.
  }
  \item{add.text}{
    A \code{logical} value controlling if the text should be added to the plot.
  }
  \item{add.contour}{
    A \code{logical} value controlling if the contour lines should be added to the plot.
  }
  \item{col}{
    A \code{character} vector of colours used in the image plot.
  }
  \item{\dots}{
    Arguments passed on.
    In \code{print.optPenaltyFusedGrid} the arguments are passed to \code{print.matrix}.
    In \code{plot.optPenaltyFusedGrid} are passed to the standard \code{plot} function.
  }
}
\value{
Invisibly returns the object (\code{x}).
}
\references{
  Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and van Wieringen, W.N. (2015).
  Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes,
  arXiv:1509.07982v1 [stat.ME].
}
\author{
  Anders Ellern Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\seealso{
\code{\link{optPenalty.fused.grid}}
}
