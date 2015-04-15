\name{optPenalty.fused}
\alias{optPenalty.fused}
\alias{optPenalty.fused.LOOCV}
\alias{optPenalty.fused.LOOCVauto}
\alias{optPenalty.fused.LOOCVgrid}

\title{
  Identify optimal ridge and fused ridge penalties
}
\description{
Functions to find the optimal ridge and fused ridge penalty parameters.
The function support leave-one-out cross-validation (LOOCV) and an approximation hereof.
Depending on the used function, general numerical optimization or a grid-based search is used.
}
\usage{
optPenalty.fused.LOOCVauto(Ylist, Tlist, lambdaF, approximate = FALSE,
                           verbose = TRUE, maxit.ridgeP.fused = 1000,
                           optimizer = "optim", maxit.optimizer = 1000,
                           debug = FALSE, ...)
optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                           lambdaMin, lambdaMax, step1 = 20,
                           lambdaFMin = lambdaMin, lambdaFMax = lambdaMax,
                           step2 = step1,
                           approximate = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{Ylist}{
   A \code{list} of matrices of data with \eqn{n_k} samples in the rows and
   \eqn{p} variables in the columns.
}
  \item{Tlist}{
   A \code{list} of the same length as \code{Ylist} of target matrices of size
   \eqn{p} times \eqn{p}.
}
  \item{lambdaF}{
   A \code{character} \code{matrix} encoding the class of penalty matrices
   to cross-validate over.
   See \code{\link{default.penalty}} for help on the construction and more
   details.
}
  \item{approximate}{
  \code{logical}. If \code{TRUE} approximate leave-one-out cross validation is
  used instead.
}
  \item{verbose}{
  \code{logical}. If \code{TRUE} additional messages are written to the terminal.
}
  \item{maxit.ridgeP.fused}{
  A \code{integer} giving the maximum number of iterations allowed for each
  fused ridge fit.
}
  \item{optimizer}{
  A \code{character} equally either \code{"optim"} or \code{"nlm"} determining
  which optimizer to use.
}
  \item{maxit.optimizer}{
  A \code{integer} giving the maximum number of iterations allowed in the
  optimization procedure.
}
  \item{debug}{
  \code{logical}. If \code{TRUE} additional output from the optimizer is appended to the
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
  Arguments passed to the optimizer in \code{optPenalty.fused.LOOCVauto} and
  arguments passed to \code{fusedRidgeS} in \code{optPenalty.fused.LOOCV}.
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
    \item \code{lambdaF} {The unique estimated optimal fused penalty matrix. I.e. the unique entries of \code{lambdaF}.}
    \item \code{lambdaF} {The estimated optimal fused penalty matrix.}
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
Ylist <- createS(n = ns, p = 6, dataset = TRUE)
Slist <- lapply(Ylist, covML)
Tlist <- default.target.fused(Slist, ns, type = "DIAES", equal = FALSE)


# Grid-based
a <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 10,
                                lambdaFMin = 1e-5, lambdaFMax = 1e3, step2 = 10,
                                approximate = FALSE, maxit = 1000)
b <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 15,
                                approximate = TRUE, maxit = 1000)



# Numerical optimization
aa <- optPenalty.fused.LOOCVauto(Ylist, Tlist, approximate = FALSE,
                                 method = "BFGS")
print(aa)
bb <- optPenalty.fused.LOOCVauto(Ylist, Tlist, approximate = TRUE,
                                 method = "BFGS")
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
Ylist <- createS(n = ns, p = 4, dataset = TRUE)
Tlist <- lapply(lapply(Ylist, covML), default.target, type = "Null")

# Construct penalty matrix
lambdaF <- default.penalty(df, type = "CartesianUnequal")

# Find optimal parameters,
# Using optim with method "Nelder-Mead"
ans1 <- optPenalty.fused.LOOCVauto(Ylist, Tlist, lambdaF = lambdaF,
                                   approximate = TRUE, verbose = FALSE)
print(ans1)
\dontrun{
# Using optim with method "BFGS"
ans2 <- optPenalty.fused.LOOCVauto(Ylist, Tlist, lambdaF = lambdaF,
                                   approximate = TRUE, verbose = FALSE,
                                   method = "BFGS")
# Using nlm
print(ans2)
ans3 <- optPenalty.fused.LOOCVauto(Ylist, Tlist, lambdaF = lambdaF,
                                   approximate = TRUE, verbose = FALSE,
                                   optimizer = "nlm")
print(ans3)
}
# Quite different results!


\dontrun{

#
# Arbitrary penalty graphs with fixed penalties!
#

# Generate some new high-dimensional data and a 2 by 2 factorial design
ns <- c(6, 5, 5, 5)
df <- expand.grid(DS = LETTERS[1:2], ER = letters[3:4])
Ylist <- createS(n = ns, p = 4, dataset = TRUE)
Tlist <- lapply(lapply(Ylist, covML), default.target, type = "Null")

lambdaF <- default.penalty(df, type = "Tensor")
print(lambdaF)  # Say we want to penalize the pair (1,2) with strength 2;
lambdaF[2,1] <- lambdaF[1,2] <- 2
print(lambdaF)

res <- optPenalty.fused.LOOCVauto(Ylist, Tlist, lambdaF = lambdaF,
                                  approximate = TRUE, optimizer = "nlm")
print(res)
}
}

