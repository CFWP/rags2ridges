\name{optFusedPenalty}
\alias{optFusedPenalty}
\alias{optFusedPenalty.LOOCVauto}
\alias{optFusedPenalty.LOOCV}

\title{
  Identify optimal ridge and fused ridge penalties
}
\description{
Functions to find the optimal ridge and fused ridge penalty parameters.
The function support leave-one-out cross-validation (LOOCV) and an approximation hereof.
Depending on the used function, general numerical optimization or a grid-based search is used.
}
\usage{
optFusedPenalty.LOOCVauto(YList, TList, Lambda,
                          approximate = FALSE, verbose = TRUE,
                          maxit.fusedRidgeS = 1000, maxit.optim = 1000,
                          optim.debug = FALSE, ...)
optFusedPenalty.LOOCV(YList, TList,
                      lambda1Min, lambda1Max, step1 = 20,
                      lambda2Min = lambda1Min, lambda2Max = lambda1Max,
                      step2 = step1,
                      approximate = FALSE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{YList}{
   A \code{list} of matrices of data with \eqn{n_k} samples in the rows and
   \eqn{p} variables in the columns.
}
  \item{TList}{
   A \code{list} of the same length as \code{YList} of target matrices of size
   \eqn{p} times \eqn{p}.
}
  \item{Lambda}{
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
  \item{maxit.fusedRidgeS}{
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
  \item{lambda1Min}{
  A single positive \code{numeric}. Minimum ridge penalty \eqn{lambda}.

}
  \item{lambda1Max}{
  A single positive \code{numeric}. Maximum ridge penalty \eqn{lambda}.
}
  \item{step1}{
  A single non-negative \code{integer}. Number grid-points used between \code{lambda1Min} and \code{lambda1Max}.
}
  \item{lambda2Min}{
  A single non-negative \code{numeric}. As for \code{lambda1Min}, but for the fused penalty assuming the complete penalty graph.
}
  \item{lambda2Max}{
  A single non-negative \code{numeric}. As for \code{lambda1Max}, but for the fused penalty assuming the complete penalty graph.
}
  \item{step2}{
  A single non-negative \code{integer}. As for \code{step1}, but for the fused penalty.
}
  \item{\dots}{
  Arguments passed to \code{optim} for \code{optFusedPenalty.LOOCVauto} and \code{fusedRidgeS} for \code{optFusedPenalty.LOOCV}.
}
}
\details{
\code{optFusedPenalty.LOOCVauto} serves a utilizes \code{\link{optim}} for identifying
the optimal fused parameters and works for general classes of penalty graphs.

\code{optFusedPenalty.LOOCV} gives a grid-based evaluation of the (approximate) LOOCV loss.
}
\value{
  \code{optFusedPenalty.LOOCVauto} returns a \code{list}:
  \itemize{
    \item \code{lambda1} {The estimated optimal ridge penalty.}
    \item \code{lambda2} {The estimated optimal fused penalty matrix if the complete penalty graph is used. If not, \code{NA} is given.}
    \item \code{LambdaP} {The estimated optimal fused penalty matrix.}
    \item \code{value} {The value of the loss function in the estimated optimum.}
  }

  \code{optFusedPenalty.LOOCV} returns a \code{list}:
  \itemize{
    \item \code{lambda1} {A \code{numeric} vector of grid values for the ridge penalty lambda1}
    \item \code{lambda2} {The \code{numeric} vector of grid values for the fused ridge penalty lambda2}
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
a <- optFusedPenalty.LOOCV(YList, TList,
                           lambda1Min = 1e-5, lambda1Max = 1e3, step1 = 10,
                           lambda2Min = 1e-5, lambda2Max = 1e3, step2 = 10,
                           approximate = FALSE, maxit = 1000)
b <- optFusedPenalty.LOOCV(YList, TList,
                           lambda1Min = 1e-5, lambda1Max = 1e3, step1 = 15,
                           approximate = TRUE, maxit = 1000)



# Numerical optimization
aa <- optFusedPenalty.LOOCVauto(YList, TList, approximate = FALSE)
print(aa)
bb <- optFusedPenalty.LOOCVauto(YList, TList, approximate = TRUE)
print(bb)

#
# Plot the results
#

# LOOCV
# Get minimums and plot
amin <- log(expand.grid(a$lambda1, a$lambda2))[which.min(a$fcvl), ]
aamin <- c(log(aa$lambda1), log(aa$lambda2))

# Plot
filled.contour(log(a$lambda1), log(a$lambda2), log(a$fcvl), color = heat.colors,
               plot.axes = {points(amin[1], amin[2], pch = 16);
                            points(aamin[1], aamin[2], pch = 16, col = "purple");
                            axis(1); axis(2)},
               xlab = "lambda1", ylab = "lambda2", main = "LOOCV")

# Approximate LOOCV
# Get minimums and plot
bmin <- log(expand.grid(b$lambda1, b$lambda2))[which.min(b$fcvl), ]
bbmin <- c(log(bb$lambda1), log(b$lambda2))

filled.contour(log(b$lambda1), log(b$lambda2), log(b$fcvl), color = heat.colors,
               plot.axes = {points(bmin[1], bmin[2], pch = 16);
                            points(bbmin[1], bbmin[2], pch = 16, col = "purple");
                            axis(1); axis(2)},
               xlab = "lambda1", ylab = "lambda2", main = "Approximate LOOCV")



#
# Arbitrary penalty graphs
#

# Generate some new high-dimensional data and a 2 by 2 factorial design
ns <- c(6, 5, 3, 2)
df <- expand.grid(Factor1 = LETTERS[1:2], Factor2 = letters[3:4])
YList <- createS(n = ns, p = 4, covariance = FALSE)
TList <- lapply(lapply(YList, covML), default.target, type = "Null")

# Construct penalty matrix
Lambda <- default.penalty(df, type = "CartesianUnequal")

# Find optimal parameters
ans <- optFusedPenalty.LOOCVauto(YList, TList,
                                 Lambda = Lambda, approximate = TRUE)
ans
}

