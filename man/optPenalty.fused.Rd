\name{optPenalty.fused}
\alias{optPenalty.fused}
\alias{optPenalty.fused.LOOCVauto}
\alias{optPenalty.fused.LOOCV}

\title{
  Identify optimal ridge and fused ridge penalties
}
\description{
Functions to find the optimal ridge and fused ridge penalty parameters.
The function support leave-one-out cross-validation (LOOCV) and an approximation hereof.
Depending on the used function, general numerical optimization or a grid-based search is used.
}
\usage{
optPenalty.fused.LOOCVauto(YList, TList, lambdaFmat,
                           approximate = FALSE, verbose = TRUE,
                           maxit.ridgeS.fused = 1000, maxit.optim = 1000,
                           optim.debug = FALSE, ...)
optPenalty.fused.LOOCV(YList, TList,
                       lambdaMin, lambdaMax, step1 = 20,
                       lambdaFMin = lambdaMin, lambdaFMax = lambdaMax,
                       step2 = step1,
                       approximate = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{YList}{
   A \code{list} of matrices of data with \eqn{n_k} samples in the rows and
   \eqn{p} variables in the columns.
}
  \item{TList}{
   A \code{list} of the same length as \code{YList} of target matrices of size
   \eqn{p} times \eqn{p}.
}
  \item{lambdaFmat}{
   A \code{character} \code{matrix} encoding the class of penalty matrices
   to cross-validate over.
   See \code{\link{default.penalty}} for help on the construction.
}
  \item{approximate}{
  \code{logical}. If \code{TRUE} approximate leave-one-out cross validation is
  used instead.
}
  \item{verbose}{
  \code{logical}. If \code{TRUE} additional messages are written to the terminal.
}
  \item{maxit.ridgeS.fused}{
  A \code{integer} giving the maximum number of iterations allowed for each
  fused ridge fit.
}
  \item{maxit.optim}{
  A \code{integer} giving the maximum number of iterations allowed in the
  optimization procedure.
}
  \item{optim.debug}{
  \code{logical}. If \code{TRUE} additional output from optim is appended to the
  output as an attribute.
}
  \item{lambdaMin}{
  A single positive \code{numeric}. Minimum ridge penalty \eqn{lambda}.

}
  \item{lambdaMax}{
  A single positive \code{numeric}. Maximum ridge penalty \eqn{lambda}.
}
  \item{step1}{
  A single non-negative \code{integer}. Number grid-points used between \code{lambdaMin} and \code{lambdaMax}.
}
  \item{lambdaFMin}{
  A single non-negative \code{numeric}. As for \code{lambdaMin}, but for the fused penalty assuming the complete penalty graph.
}
  \item{lambdaFMax}{
  A single non-negative \code{numeric}. As for \code{lambdaMax}, but for the fused penalty assuming the complete penalty graph.
}
  \item{step2}{
  A single non-negative \code{integer}. As for \code{step1}, but for the fused penalty.
}
  \item{\dots}{
  Arguments passed to \code{optim} for \code{optPenalty.fused.LOOCVauto} and \code{fusedRidgeS} for \code{optPenalty.fused.LOOCV}.
}
}
\details{
\code{optPenalty.fused.LOOCVauto} serves a utilizes \code{\link{optim}} for identifying
the optimal fused parameters and works for general classes of penalty graphs.

\code{optPenalty.fused.LOOCV} gives a grid-based evaluation of the (approximate) LOOCV loss.
}
\value{
  \code{optPenalty.fused.LOOCVauto} returns a \code{list}:
  \itemize{
    \item \code{lambda} {The estimated optimal ridge penalty.}
    \item \code{lambdaF} {The estimated optimal fused penalty matrix if the complete penalty graph is used. If not, \code{NA} is given.}
    \item \code{lambdaFmat} {The estimated optimal fused penalty matrix.}
    \item \code{value} {The value of the loss function in the estimated optimum.}
  }

  \code{optPenalty.fused.LOOCV} returns a \code{list}:
  \itemize{
    \item \code{lambda} {A \code{numeric} vector of grid values for the ridge penalty \code{lambda}}
    \item \code{lambdaF} {The \code{numeric} vector of grid values for the fused ridge penalty \code{lambdaF}}
    \item \code{fcvl} {The \code{numeric} \code{matrix} of evaluations of the loss function}
  }
}
\references{
Bilgrau, AE; Peeters CFW; Eriksen, PS; Boegsted, M; & van Wieringen, WN
(in preparation) "Fused Ridge Estimation of Multiple Inverse Covariance
Matrices from High-Dimensional Data Classes"
}
\author{
Anders Ellern Bilgrau \cr
Carel F.W. Peeters <cf.peeters@vumc.nl> \cr
Wessel N. van Wieringen
}
\seealso{
See also \code{\link{default.penalty}}, \code{optPenalty.LOOCV}.
}
\examples{
# Generate some (not so) high-dimensional data witn (not so) many samples
ns <- c(4, 5, 6)
YList <- createS(n = ns, p = 6, covariance = FALSE)
SList <- lapply(YList, covML)
TList <- lapply(SList, default.target, type = "DIAES")


# Grid-based
a <- optPenalty.fused.LOOCV(YList, TList,
                            lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 10,
                            lambdaFMin = 1e-5, lambdaFMax = 1e3, step2 = 10,
                            approximate = FALSE, maxit = 1000)
b <- optPenalty.fused.LOOCV(YList, TList,
                            lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 15,
                            approximate = TRUE, maxit = 1000)



# Numerical optimization
aa <- optPenalty.fused.LOOCVauto(YList, TList, approximate = FALSE)
print(aa)
bb <- optPenalty.fused.LOOCVauto(YList, TList, approximate = TRUE)
print(bb)

#
# Plot the results
#

# LOOCV
# Get minimums and plot
amin <- log(expand.grid(a$lambda, a$lambdaF))[which.min(a$fcvl), ]
aamin <- c(log(aa$lambda), log(aa$lambdaF))

# Plot
filled.contour(log(a$lambda), log(a$lambdaF), log(a$fcvl), color = heat.colors,
               plot.axes = {points(amin[1], amin[2], pch = 16);
                            points(aamin[1], aamin[2], pch = 16, col = "purple");
                            axis(1); axis(2)},
               xlab = "lambda", ylab = "lambdaF", main = "LOOCV")

# Approximate LOOCV
# Get minimums and plot
bmin <- log(expand.grid(b$lambda, b$lambdaF))[which.min(b$fcvl), ]
bbmin <- c(log(bb$lambda), log(b$lambdaF))

filled.contour(log(b$lambda), log(b$lambdaF), log(b$fcvl), color = heat.colors,
               plot.axes = {points(bmin[1], bmin[2], pch = 16);
                            points(bbmin[1], bbmin[2], pch = 16, col = "purple");
                            axis(1); axis(2)},
               xlab = "lambda", ylab = "lambdaF", main = "Approximate LOOCV")



#
# Arbitrary penalty graphs
#

# Generate some new high-dimensional data and a 2 by 2 factorial design
ns <- c(6, 5, 3, 2)
df <- expand.grid(Factor1 = LETTERS[1:2], Factor2 = letters[3:4])
YList <- createS(n = ns, p = 4, covariance = FALSE)
TList <- lapply(lapply(YList, covML), default.target, type = "Null")

# Construct penalty matrix
lambdaFmat <- default.penalty(df, type = "CartesianUnequal")

# Find optimal parameters
ans <- optPenalty.fused.LOOCVauto(YList, TList,
                                  lambdaFmat = lambdaFmat, approximate = TRUE)
ans
}

