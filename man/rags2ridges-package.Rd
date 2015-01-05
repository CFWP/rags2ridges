\name{rags2ridges-package}
\alias{rags2ridges-package}
\alias{rags2ridges}
\docType{package}
\title{
Ridge estimation for high-dimensional precision matrices
}
\description{
Package contains proper L2-penalized ML estimators for the precision matrix as well as supporting functions to employ these estimators in a graphical modeling setting.
}
\details{
The main function of the package is \code{\link{ridgeS}} which enables archetypal and proper alternative ML ridge estimation of the precision matrix. The alternative 
ridge estimators can be found in van Wieringen and Peeters (2014) and encapsulate both target and non-target shrinkage for the multivariate normal precision 
matrix. The estimators are analytic and enable estimation in large \eqn{p} small \eqn{n} settings. Supporting functions to employ these estimators in a graphical
modeling setting are also given. These supporting functions enable, a.o., the determination of the optimal value of the penalty parameter, the determination of 
the support of a shrunken precision estimate, as well as various visualization options.

Function overview:
\itemize{
\item Function for (proper) ridge estimation of the precision matrix
	\itemize{
	\item \code{\link{ridgeS}}
	}
\item Functions for penalty parameter selection
	\itemize{
	\item \code{\link{conditionNumberPlot}}
	\item \code{\link{optPenalty.aLOOCV}}
	\item \code{\link{optPenalty.LOOCV}}
	\item \code{\link{optPenalty.LOOCVauto}}
	}
\item Functions for loss/entropy/fit evaluation
	\itemize{
	\item \code{\link{evaluateSfit}}
	\item \code{\link{KLdiv}}
	\item \code{\link{loss}}
	}
\item Functions for block-independence testing
	\itemize{
	\item \code{\link{GGMblockNullPenalty}}
	\item \code{\link{GGMblockTest}}
	}
\item Function for support determination
	\itemize{
	\item \code{\link{sparsify}}
	}
\item Functions for (network) visualization
	\itemize{
	\item \code{\link{edgeHeat}}
	\item \code{\link{ridgePathS}}
	\item \code{\link{Ugraph}}
	}
\item Functions for topology statistics
	\itemize{
	\item \code{\link{GGMnetworkStats}}
	\item \code{\link{GGMpathStats}}
	}
\item Wrapper function
	\itemize{
	\item \code{\link{fullMontyS}}
	}
\item Support functions
	\itemize{
	\item \code{\link{adjacentMat}}
	\item \code{\link{covML}}
	\item \code{\link{default.target}}
	\item \code{\link{evaluateS}}
	\item \code{\link{pcor}}
	\item \code{\link{symm}}
	}
}
}
\author{
Carel F.W. Peeters, Wessel, N. van Wieringen
\cr
Maintainer: Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\references{
van Wieringen, W.N. and Peeters, C.F.W. (2014). Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data. arXiv:1403.0904 [stat.ME].
}
