\name{rags2ridges-package}
\alias{rags2ridges-package}
\alias{rags2ridges}
\docType{package}
\title{
Ridge estimation for high-dimensional precision matrices
}
\description{
Package contains proper L2-penalized ML estimators for the precision matrix as well as supporting functions to employ these estimators in a (integrative or meta-analytic) graphical modeling setting.
}
\details{
The main function of the package is \code{\link{ridgeP}} which enables archetypal and proper alternative ML ridge estimation of the precision matrix. The alternative
ridge estimators can be found in van Wieringen and Peeters (2015) and encapsulate both target and non-target shrinkage for the multivariate normal precision
matrix. The estimators are analytic and enable estimation in large \eqn{p} small \eqn{n} settings. Supporting functions to employ these estimators in a graphical
modeling setting are also given. These supporting functions enable, a.o., the determination of the optimal value of the penalty parameter, the determination of
the support of a shrunken precision estimate, as well as various visualization options.

The package has a modular setup. The \emph{core module} (rags2ridges.R) contains the functionality stated above.
The \emph{fused module} (rags2ridgesFused.R) extends the functionality of the core module to the joint estimation of
multiple precision matrices from (aggregated) high-dimensional data consisting of distinct classes.
The result is a targeted fused ridge estimator that is of use when the precision matrices of the
constituent classes are believed to chiefly share the same structure while potentially differing
in a number of locations of interest. The fused module also contains supporting functions for
integrative or meta-analytic Gaussian graphical modeling. The third module is the
\emph{miscellaneous module} (rags2RidgesMisc.R) which contains assorted hidden functions.

Function overview \emph{core module}:
\itemize{
\item Function for (proper) ridge estimation of the precision matrix
	\itemize{
	\item \code{\link{ridgeP}}
	}
\item Functions for penalty parameter selection
	\itemize{
	\item \code{\link{CNplot}}
	\item \code{\link{optPenalty.aLOOCV}}
	\item \code{\link{optPenalty.kCV}}
	\item \code{\link{optPenalty.kCVauto}}
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
	\item \code{\link{GGMmutualInfo}}
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
	\item \code{\link{covMLknown}}
	\item \code{\link{default.target}}
	\item \code{\link{evaluateS}}
	\item \code{\link{pcor}}
	\item \code{\link{symm}}
	}
}

Function overview \emph{fused module}:
\itemize{
\item Function for targeted fused ridge estimation of multiple precision matrices
	\itemize{
	\item \code{\link{ridgeP.fused}}
	}
\item Function for fused penalty parameter selection
	\itemize{
	\item \code{\link{optPenalty.fused}}
	}
\item Functions for loss/entropy/fit evaluation
	\itemize{
	\item \code{\link{KLdiv.fused}}
	\item \code{\link{NLL}}
	}
\item Function for testing the necessity of fusion
	\itemize{
	\item \code{\link{fused.test}}
	}
\item Function for support determination
	\itemize{
	\item \code{\link{sparsify.fused}}
	}
\item Functions for topology statistics
	\itemize{
	\item \code{\link{GGMnetworkStats.fused}}
	\item \code{\link{GGMpathStats.fused}}
	}
\item Support functions
	\itemize{
	\item \code{\link{createS}}
	\item \code{\link{default.penalty}}
	\item \code{\link{default.target.fused}}
	\item \code{\link{getKEGGPathway}}
	\item \code{\link{isSymmetricPD}}
	\item \code{\link{is.Xlist}}
	\item \code{\link{kegg.target}}
	\item \code{\link{plot.ptest}}
	\item \code{\link{pooledS}}
	\item \code{\link{print.optPenaltyFusedGrid}}
	\item \code{\link{print.ptest}}
	\item \code{\link{rmvnormal}}
	}
}

Calls of interest to \emph{miscellaneous module}:
\itemize{
  \item \code{rags2ridges:::.TwoCents()} ~~(Unsolicited advice)
  \item \code{rags2ridges:::.Brooke()} ~~(Endorsement)
  \item \code{rags2ridges:::.JayZScore()} ~~(The truth)
  \item \code{rags2ridges:::.theHoff()} ~~(Wish)
  \item \code{rags2ridges:::.rags2logo()} ~~(Warm welcome)
}
}
\author{
Carel F.W. Peeters, Anders Ellern Bilgrau, Wessel, N. van Wieringen
\cr
Maintainer: Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\references{
  Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and van Wieringen, W.N. (2015).
  Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes,
  arXiv:1509.07982v1 [stat.ME].

  Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2016).
  The Spectral Condition Number Plot for Regularization Parameter Determination,
  arXiv:1608.04123v1 [stat.CO].

  van Wieringen, W.N. & Peeters, C.F.W. (2016).
  Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data,
  Computational Statistics & Data Analysis, vol. 103: 284-303.
  Also available as arXiv:1403.0904v3 [stat.ME].

  van Wieringen, W.N. & Peeters, C.F.W. (2015).
  Application of a New Ridge Estimator of the Inverse Covariance Matrix
  to the Reconstruction of Gene-Gene Interaction Networks.
  In: di Serio, C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.)
  `Computational Intelligence Methods for Bioinformatics and Biostatistics'.
  Lecture Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
}
