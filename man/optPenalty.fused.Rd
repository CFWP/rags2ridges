\name{optPenalty.fused}
\alias{optPenalty.fused}
\alias{optPenalty.fused.LOOCV}
\alias{optPenalty.fused.LOOCVauto}
\alias{optPenalty.fused.LOOCVgrid}

\title{
Identify optimal ridge and fused ridge penalties
}
\description{
Functions to find the optimal ridge and fusion penalty parameters via
leave-one-out cross validation.
The functions support leave-one-out cross-validation (LOOCV), \eqn{k}-fold CV,
and two forms of approximate LOOCV.
Depending on the used function, general numerical optimization or a grid-based
search is used.
}
\usage{
optPenalty.fused(Ylist, Tlist,
                 lambda = default.penalty(Ylist),
                 cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
                 k = 10, grid = FALSE, ...)

# A simple wrapper for:
optPenalty.fused.LOOCVauto(Ylist, Tlist, lambda,
                           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
                           k = 10, verbose = TRUE, lambda.init,
                           maxit.ridgeP.fused = 1000,
                           optimizer = "optim", maxit.optimizer = 1000,
                           debug = FALSE,
                           optim.control = list(trace = verbose, maxit = maxit.optimizer),
                           ...)
optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                           lambdaMin, lambdaMax, step1 = 20,
                           lambdaFMin = lambdaMin, lambdaFMax = lambdaMax,
                           step2 = step1,
                           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
                           k = 10,
                           verbose = TRUE, ...)
}
\arguments{
  \item{Ylist}{
    A \code{list} of \eqn{G} matrices of data with \eqn{n_g} samples in the rows
    and \eqn{p} variables in the columns corresponding to \eqn{G} classes of
    data.
  }
  \item{Tlist}{
    A \code{list} of \eqn{G} of p.d. class target matrices of size
    \eqn{p} times \eqn{p}.
  }
  \item{lambda}{
    A symmetric \code{character} \code{matrix} encoding the class of penalty
    matrices to cross-validate over.
    The diagonal elements correspond to the class-specific ridge penalties
    whereas the off-diagonal elements correspond to the fusion penalties.
    The unique elements of lambda specify the penalties to determine by
    the method specified by \code{cv.method}.
    The penalties can be fixed if they are coercible to numeric values, such as
    e.g. \code{"0"}, \code{"2.71"} or \code{"3.14"}.
    Fusion between pairs can be "left out"" using either of \code{""},
    \code{NA}, \code{"NA"}, or \code{"0"}.
    See \code{\link{default.penalty}} for help on the construction hereof and
    more details.
    Unused and can be omitted if \code{grid == TRUE}.
  }
  \item{cv.method}{
    \code{character} giving the cross-validation (CV) to use. The allowed values
    are \code{"LOOCV"}, \code{"aLOOCV"}, \code{"sLOOCV"}, \code{"kCV"} for
    leave-one-out cross validation (LOOCV), appproximate LOOCV, special LOOCV,
    and k-fold CV, respectively.
  }
  \item{k}{
    \code{integer} giving the number of approximately equally sized parts each
    class is partioned into for \eqn{k}-fold CV.
    Only use if \code{cv.method} is \code{"kCV"}.
  }
  \item{verbose}{
    \code{logical}. If \code{TRUE}, progress information is printed to the
    console.
  }
  \item{lambda.init}{
    A \code{numeric} penalty \code{matrix} of initial values passed to the
    optimizer. If omitted, the function selects a starting values using
    a common ridge penaltiy (determined by 1D optimization)
    and all fusion penalties to zero.
  }
  \item{maxit.ridgeP.fused}{
    A \code{integer} giving the maximum number of iterations allowed for each
    fused ridge fit.
  }
  \item{optimizer}{
    \code{character}. Either \code{"optim"} or \code{"nlm"} determining
    which optimizer to use.
  }
  \item{maxit.optimizer}{
    A \code{integer} giving the maximum number of iterations allowed in the
    optimization procedure.
  }
  \item{debug}{
    \code{logical}. If \code{TRUE} additional output from the optimizer is
    appended to the output as an attribute.
  }
  \item{lambdaMin}{
    A single positive \code{numeric}. Minimum ridge penalty \eqn{lambda}.
  }
  \item{lambdaMax}{
    A single positive \code{numeric}. Maximum ridge penalty \eqn{lambda}.
  }
  \item{step1}{
    A single non-negative \code{integer}. Number grid-points used between
    \code{lambdaMin} and \code{lambdaMax}.
  }
  \item{lambdaFMin}{
    A single non-negative \code{numeric}. As for \code{lambdaMin}, but for the
    fused penalty assuming the complete penalty graph.
  }
  \item{lambdaFMax}{
    A single non-negative \code{numeric}. As for \code{lambdaMax}, but for the
    fused penalty assuming the complete penalty graph.
  }
  \item{step2}{
    A single non-negative \code{integer}. As for \code{step1}, but for the fused
    penalty.
  }
  \item{grid}{
    \code{logical.} Should a grid based search be used? Default is \code{FALSE}.
  }
  \item{optim.control}{
    A \code{list} of control arguments for \code{\link{optim}}.
  }
  \item{\dots}{
    For \code{optPenalty.fused}, arguments are passed to
    \code{optPenalty.fused.LOOCVgrid} or \code{optPenalty.fused.LOOCVauto}
    depending on the value of \code{grid}.
    In \code{optPenalty.fused.LOOCVgrid}, arguments are passed to
    \code{fusedRidgeS}.
    In \code{optPenalty.fused.LOOCVauto}, arguments are passed to the optimizer.
  }
}
\details{
\code{optPenalty.fused.LOOCVauto} serves a utilizes \code{\link{optim}} for
identifying the optimal fused parameters and works for general classes of
penalty graphs.

\code{optPenalty.fused.LOOCV} gives a grid-based evaluation of the (approximate)
LOOCV loss.
}
\value{
\code{optPenalty.fused.LOOCVauto} returns a \code{list}:\cr
\item{Plist}{A \code{list} of the precision estimates for the optimal
  parameters.}
\item{lambda}{The estimated optimal fused penalty matrix.}
\item{lambda.unique}{The unique entries of the \code{lambda}.
  A more concise overview of \code{lambda}}
\item{value}{The value of the loss function in the estimated optimum.}

\code{optPenalty.fused.LOOCV} returns a \code{list}:\cr
\item{ridge}{A \code{numeric} vector of grid values for the ridge
  penalty}
\item{fusion}{The \code{numeric} vector of grid values for the fusion
  penalty}
\item{fcvl}{The \code{numeric} \code{matrix} of evaluations of the
  loss function}
}
\references{
Bilgrau, AE; Peeters CFW; Eriksen, PS; Boegsted, M; & van Wieringen, WN
(In preparation) "Targeted Fused Ridge Estimation of Multiple Inverse Covariance
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
Tlist <- default.target.fused(Slist, ns, type = "DIAES")


# Grid-based
a <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 10,
                                lambdaFMin = 1e-5, lambdaFMax = 1e3, step2 = 10,
                                cv.method = "LOOCV", maxit = 1000)
b <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 15,
                                cv.method = "aLOOCV", maxit = 1000)
c <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 15,
                                cv.method = "sLOOCV", maxit = 1000)
d <- optPenalty.fused.LOOCVgrid(Ylist, Tlist,
                                lambdaMin = 1e-5, lambdaMax = 1e3, step1 = 15,
                                cv.method = "kCV", k = 2, maxit = 1000)

# Numerical optimization (uses the default "optim" optimizer with method "BFGS")
aa <- optPenalty.fused.LOOCVauto(Ylist, Tlist, cv.method = "LOOCV",
                                 method = "BFGS")
print(aa)
bb <- optPenalty.fused.LOOCVauto(Ylist, Tlist, cv.method = "aLOOCV",
                                 method = "BFGS")
print(bb)
cc <- optPenalty.fused.LOOCVauto(Ylist, Tlist, cv.method = "sLOOCV",
                                 method = "BFGS")
print(cc)
dd <- optPenalty.fused.LOOCVauto(Ylist, Tlist, cv.method = "kCV", k = 3,
                                 method = "BFGS")
print(dd)


#
# Plot the results
#

# LOOCV
# Get minimums and plot
amin <- log(expand.grid(a$lambda, a$lambdaF))[which.min(a$fcvl), ]
aamin <- c(log(aa$lambda[1,1]), log(aa$lambda[1,2]))

# Plot
filled.contour(log(a$lambda), log(a$lambdaF), log(a$fcvl), color = heat.colors,
               plot.axes = {points(amin[1], amin[2], pch = 16);
                            points(aamin[1], aamin[2], pch = 16, col = "purple");
                            axis(1); axis(2)},
               xlab = "lambda", ylab = "lambdaF", main = "LOOCV")

# Approximate LOOCV
# Get minimums and plot
bmin <- log(expand.grid(b$lambda, b$lambdaF))[which.min(b$fcvl), ]
bbmin <- c(log(bb$lambda[1,1]), log(unique(bb$lambda[1,2])))

filled.contour(log(b$lambda), log(b$lambdaF), log(b$fcvl), color = heat.colors,
               plot.axes = {points(bmin[1], bmin[2], pch = 16);
                            points(bbmin[1], bbmin[2], pch = 16, col ="purple");
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
lambda <- default.penalty(df, type = "CartesianUnequal")

# Find optimal parameters,
# Using optim with method "Nelder-Mead" with "special" LOOCV
ans1 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
                         cv.method = "sLOOCV", verbose = FALSE)
print(ans1$lambda.unique)

# By approximate LOOCV using optim with method "BFGS"
ans2 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
                         cv.method = "aLOOCV", verbose = FALSE,
                         method = "BFGS")
print(ans2$lambda.unique)

# By LOOCV using nlm
ans3 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
                         cv.method = "LOOCV", verbose = FALSE,
                         optimizer = "nlm")
print(ans3$lambda.unique)

# Quite different results!


#
# Arbitrary penalty graphs with fixed penalties!
#

# Generate some new high-dimensional data and a 2 by 2 factorial design
ns <- c(6, 5, 5, 5)
df <- expand.grid(DS = LETTERS[1:2], ER = letters[3:4])
Ylist <- createS(n = ns, p = 4, dataset = TRUE)
Tlist <- lapply(lapply(Ylist, covML), default.target, type = "Null")

lambda <- default.penalty(df, type = "Tensor")
print(lambda)  # Say we want to penalize the pair (1,2) with strength 2.1;
lambda[2,1] <- lambda[1,2] <- 2.1
print(lambda)

# Specifiying starting values is also possible:
init <- diag(length(ns))
init[2,1] <- init[1,2] <- 2.1

res <- optPenalty.fused(Ylist, Tlist, lambda = lambda, lambda.init = init,
                        cv.method = "aLOOCV", optimizer = "nlm")
print(res)
}

