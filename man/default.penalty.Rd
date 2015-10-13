\name{default.penalty}
\alias{default.penalty}
\title{
  Construct commonly used penalty matrices
}
\description{
  Function that constructs default or commonly use penalty matrices
  according to a (factorial) study design.
  The constructed penalty matrix can be used directly in
  \code{\link{optPenalty.fused.auto}} or serve as basis for modification.
}
\usage{
  default.penalty(G, df,
                  type = c("Complete", "CartesianEqual",
                           "CartesianUnequal", "TensorProd"))
}
\arguments{
  \item{G}{
    A \code{numeric} giving the number of classes. Can also be a \code{list} of
    length \code{G} such as the usual argument \code{Slist} from other
    \pkg{rags2ridges} functions.
    Can be omitted if \code{df} is given.
  }
  \item{df}{
    A \code{data.frame} with \code{G} rows and factors in the columns.
    Note, the columns has to be of type \code{factor}.
    Can be omitted when \code{G} is given and \code{type == "Complete"}.
    The factors can be ordered.
  }
  \item{type}{
    A character giving the type of fused penalty graph to construct.
    Should be \code{'Complete'} (default), \code{'CartesianEqual'}, or
    \code{'CartesianUnequal'} or \code{'TensorProd'} or an unique abbreviation
   hereof. See details.
  }
}
\details{
 The \code{type} gives a number of common choices for the penalty matrix:
 \itemize{
   \item \code{'Complete'} is the complete penalty graph with equal
     penalties.
   \item \code{'CartesianEqual'} corresponds to a penalizing along each
     "direction" of factors with a common penalty. The choice is named
     Cartesian as it is the Cartesian graph product of the complete penalty
     graphs for the individual factors.
   \item \code{'CartesianUnequal'} corresponds to a penalizing each
     direction of factors with individual penalties.
   \item \code{'TensorProd'} correspond to penalizing the "diagonals" only.
     It is equivalent to the graph tensor products of the complete graphs
     for each individual factor.
  }
}
\value{
  Returns a \code{G} by \code{G} character matrix which specify the class
  of penalty graphs to be used.
  The output is suitable as input for the penalty matrix
  used in \code{\link{optPenalty.fused.auto}}.
}
\references{
  Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and van Wieringen, W.N. (2015).
  Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes,
  arXiv:1509.07982v1 [stat.ME].
}
\author{
  Anders E. Bilgrau,
  Carel F.W. Peeters <cf.peeters@vumc.nl>,
  Wessel N. van Wieringen
}
\seealso{
  \code{\link{ridgeP.fused}},
  \code{\link{optPenalty.fused}},
  \code{\link{default.target}}
}
\examples{
  # Handling one-way designs
  default.penalty(2)
  default.penalty(4)
  Slist <- vector("list", 6)
  default.penalty(Slist)   # The function uses only the length of the list
  df0 <- expand.grid(Factor = c("lvl1", "lvl2"))
  default.penalty(df0)

  # A more elaborate example
  df1 <- expand.grid(DS = c("DS1", "DS2", "DS3"), ER = c("ER+", "ER-"))

  # Usage (various interface demonstrations)
  default.penalty(6, df1, type = "Complete")
  default.penalty(6, type = "CartesianEqual")  # GIVES WARNING
  default.penalty(6, df1, type = "CartesianEqual")
  default.penalty(Slist, df1, type = "CartesianEqual")
  default.penalty(6, df1, type = "CartesianUnequal")
  default.penalty(df1)

  # A 2 by 2 by 2 design
  df2 <- expand.grid(A = c("A1", "A2"), B = c("B1", "B2"), C = c("C1", "C3"))
  default.penalty(df2)
  default.penalty(df2, type = "CartesianEqual")
  default.penalty(df2, type = "CartesianUnequal")
}


