\name{NEWS}
\title{News for Package \pkg{rags2ridges}}
\newcommand{\cpkg}{\href{http://cran.r-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 2.0}{
  \itemize{
    \item Added this \code{NEWS} file!
    \item Added the fused ridge module. This provides estimation functionality
      for estimation of multiple precision matrices for multiple classes.
      Functions from this module are generally suffixed with \code{.fused}.
    \item Argument \code{verticle} in various functions have been renamed to
      \code{vertical}. Sorry for any inconvinience.
    \item Added LOOCV and approximations hereof to determine optimal ridge and
      fused penalties. The fused pairs to cross-validate over are specified by
      a character matrix.
    \item \cpkg{rags2ridges} now uses \cpkg{Rcpp} and \cpkg{RcppArmadillo} with
      core functions written in \code{C++}. The package should now be at least
      two orders of magnitude faster in most cases.
  }
}

\section{Changes in version 1.4}{
  \itemize{
    \item Inclusion hidden function \code{.pathContribution} for usage in
      \code{GGMpathStats} function
    \item Inclusion hidden function \code{.path2string} for usage in
      \code{GGMpathStats} function
    \item Inclusion hidden function \code{.pathAndStats} for usage in
      \code{GGMpathStats} function
    \item Inclusion hidden function \code{.cvl} for usage in
      \code{optPenalty.LOOCVauto} function
    \item Inclusion hidden function \code{.lambdaNullDist} for usage in
      \code{GGMblockNullPenalty} function
    \item Inclusion hidden function \code{.blockTestStat} for usage in
      \code{GGMblockTest} function
    \item Inclusion function that expresses the covariance between a pair of
      variables as a sum of path weights: \code{GGMpathStats}
    \item Inclusion function that determines the optimal penalty parameter
      value by application of the Brent algorithm to the LOOCV log-likelihood:
      \code{optPenalty.LOOCVauto}
    \item Inclusion function that generates the distribution of the penalty
      parameter under the null hypothesis of block independence:
      \code{GGMblockNullPenalty}
    \item Inclusion function that performs a permutation test for block
      structure in the precision matrix: \code{GGMblockTest}
    \item Inclusion wrapper function: \code{fullMontyS}
    \item New features updated \code{optPenalty.aLOOCV} function:
      For scalar matrix targets the complete solution path depends on only
      1 eigendecomposition and 1 matrix inversion.
      Meaning: the function is sped up somewhat by lifting redundant
      inversions out of \code{for} loops
      Optional graph now plots the approximated LOOCV negative log-likelihood
      instead of ln(approximated LOOCV negative log-likelihood).
      Legend in optional graph has been adapated accordingly
    \item New features updated \code{optPenalty.LOOCV} function:
      Optional graph now plots the LOOCV negative log-likelihood instead of
      ln(LOOCV negative log-likelihood).
      Legend in optional graph has been adapated accordingly
    \item New features updated \code{default.target} function:
      Inclusion new default target option: \code{type = DIAES}. Gives diagonal
      matrix with inverse of average of eigenvalues of S as entries
    \item New features updated \code{GGMnetworkStats} function:
      Now also assesses (and returns a logical) if graph/network is chordal
      Now also includes assesment of the eigenvalue centrality
      Now includes option to have list or table output
    \item New features updated \code{ridgePathS} function:
      Sped up considerably for rotation equivariant alternative estimator
      By avoidance of redundant eigendecompositions and inversions
      Now catches breakdown due to rounding preculiarities when
      \code{plotType = "pcor"}
    \item New features updated \code{sparsify} function:
      Inclusion new thresholding function \code{top}: retainment of top elements
      based on absolute partial correlation
      Inclusion output option: When \code{output = "light"}, only the (matrix)
      positions of the zero and non-zero elements are returned
      Function no longer dependent on GeneNet; now makes direct use of
      \cpkg{fdrtool}.
      Function now also prints some general information on the number of edges
      retained.
    \item Corrected small error in \code{evaluateSfit} function
      The \code{dir} argument was not properly used previously.
  }
}

\section{Changes in version 1.3}{
  \itemize{
    \item Inclusion hidden function \code{.ridgeSi} for usage in
      \code{conditionNumberPlot} function.
    \item Inclusion hidden function \code{.eigShrink} for usage in (a.o.)
      \code{ridgeS} function.
    \item Inclusion function calculating various network statistics from a
      sparse matrix: \code{GGMnetworkStats}
    \item Inclusion function for visual inspection fit of regularized precision
      matrix to sample covariance matrix: \code{evaluateSfit}
    \item Inclusion function for visualization of regularization paths:
      \code{ridgePathS}
    \item Inclusion function for default target matrix generation:
      \code{default.target}
    \item New features updated \code{evaluateS} function:
      The printed output of the \code{evaluateS} function is now aligned
      Calculation spectral condition number has been improved
    \item \code{conditionNumber} function now called \code{conditionNumberPlot}.
      Updated function has new features:
      Main plot can now be obtained with either the spectral (l2) or the
      (approximation to) l1 condition number
      Main plot can now be amended with plot of the relative distance to the
      set of singular matrices
      The title of the main plot can now be suppressed
      One can now obtain numeric output from the function: lambdas and
      condition numbers
    \item New features updated \code{sparsify} function:
      Changed \code{type = c("threshold", "localFDR")} to
      \code{threshold = c("absValue", "localFDR")} (clarifying nomenclature)
      Changed \code{threshold} to \code{absValueCut} (clarifying nomenclature)
      Will now output both sparsified partial correlation/standardized
      precision matrix and the sparsified precison matrix >>
      when input consists of the unstandardized precision matrix
    \item New features updated \code{ridgeS} function:
      Contains an improved evaluation of the target matrix possibly being a
      null matrix
      Now evaluates if a rotation equivariant alternative estimator ensues for
      a given target matrix
      When rotation equivariant alternative estimator ensues, computation is
      sped up considerably by circumventing the matrix square root
    \item \code{optPenaltyCV} function now called \code{optPenalty.LOOCV}, for sake
      of (naming) consistency. Updated function has new features:
      \code{targetScale} option has been removed
      Replaced \code{log} in optional graph by \code{ln}
      Visual layout of optional graph now more in line with recommendations by
      Tufte.
    \item New features updated \code{optPenalty.aLOOCV} function:
      Replaced \code{log} in optional graph by \code{ln}
      Visual layout of optional graph now more in line with recommendations by Tufte
    \item Computation optimal penalty in \code{conditionNumberPlot},
      \code{optPenalty.aLOOCV} and \code{optPenalty.LOOCV} functions sped up
      considerably: For rotation equivariant alternative estimator
      By usage new ridgeS and avoidance of redundant eigendecompositions
    \item Default target in \code{ridgeS}, \code{conditionNumberPlot},
      \code{optPenalty.aLOOCV} and \code{optPenalty.LOOCV} now \code{"DAIE"}
      option in \code{default.target}.
  }
}


\section{Changes in version 1.2}{
  \itemize{
    \item Inclusion function for ML estimation of the sample covariance matrix:
      \code{covML}
    \item Inclusion function for approximate leave-one-out cross-validation:
      \code{optPenalty.aLOOCV}
    \item Inclusion function \code{conditionNumber} to visualize the
      spectral condition number over the regularization path
    \item Inclusion function \code{evaluateS} to evaluate basic
      properties of a covariance matrix
    \item Inclusion function \code{KLdiv} that calculates the
      Kullback-Leibler divergence between two normal distributions
    \item Inclusion option to suppress on-screen output in \code{sparsify}
      function
    \item Corrected small error in \code{optPenaltyCV} function
    \item Both \code{optPenaltyCV} and \code{optPenalty.aLOOCV} now utilize
      \code{covML} instead of \code{cov}
    \item Default output option in \code{optPenaltyCV} (as in
      \code{optPenalty.aLOOCV}) is now \code{light}
  }
}





