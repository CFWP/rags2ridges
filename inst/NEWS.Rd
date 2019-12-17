\name{NEWS}
\newcommand{\cpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}}

\title{News for Package \cpkg{rags2ridges}}

\section{Changes in version 2.2.2}{
  \subsection{Notifications:}{
    \itemize{
      \item Hot fix due to class changes in "matrix". No major visible user changes.
      \item \code{CNplot} function again updated: higher max. iterations for Lanczos method
    }
  }
}

\section{Changes in version 2.2.1}{
  \subsection{Notifications:}{
    \itemize{
      \item Hot fix due to new RNG. No major visible user changes.
      \item \code{CNplot} function updated
    }
  }
}

\section{Changes in version 2.2}{
  \subsection{Notifications:}{
    \itemize{
      \item \code{optPenalty.LOOCV} is deprecated. Please use \code{optPenalty.kCV} instead
      \item \code{optPenalty.LOOCVauto} is deprecated. Please use \code{optPenalty.kCVauto} instead
    }
  }
}


\section{Changes in version 2.1.1}{
  \subsection{Documentation:}{
    \itemize{
      \item Updated \code{CITATION} file
      \item Updated \code{README} file
    }
  }
  \subsection{Adjustments:}{
    \itemize{
      \item \code{sparsify} now has an additional thresholding option: 'connected'
    }
  }
}


\section{Changes in version 2.1}{
  \subsection{Documentation:}{
    \itemize{
      \item Updated \code{CITATION} file
      \item Updated \code{README} file
    }
  }
  \subsection{Bug fixes:}{
    \itemize{
      \item Fixed bug in \code{Ugraph}:
        \itemize{
          \item Incorrectly stated before that all igraph layouts were supported.
          \item Now they indeed are supported.
      }
    }
  }
  \subsection{Notifications:}{
    \itemize{
      \item \code{conditionNumberPlot} is deprecated. Please use \code{CNplot} instead
      \item Features of the \code{CNplot} function (above and beyond \code{conditionNumberPlot}):
        \itemize{
          \item The \code{digitLoss} and \code{rlDist} arguments have been removed
          \item These arguments have been replaced with the logical argument \code{Iaids}
          \item \code{Iaids = TRUE} amends the basic condition number plot with interpretational aids
          \item These aids are the approximate loss in digits of accuracy and and approximation of the acceleration             along the regularization path of the condition number
          \item Argument \code{main} is now a character argument
          \item Argument \code{value} now by default takes the value 1e-100 (convenient)
          \item Now uses C++ functionalty for additional speed
      }
    }
  }
  \subsection{Adjustments:}{
    \itemize{
      \item \code{edgeHeat} now has aligned x-axis labels
      \item The visualizations of the \code{optPenalty.LOOCV} and \code{optPenalty.aLOOCV} functions
      now will no longer produce horizontal and/or vertical lines that fall outside the boundaries
      of the figure
      \item \code{optPenalty.LOOCV} now uses log-equidistant penalty grid for optimal penalty
      parameter determination (this also enhances the visualization)
      \item New features updated \code{optPenalty.aLOOCV} function:
        \itemize{
          \item Function has been sped up by killing redundant inversion
          \item now uses log-equidistant penalty grid for optimal penalty
          parameter determination (this also enhances the visualization)
      }
      \item New features updated \code{Ugraph} function:
        \itemize{
          \item One can now also specify vertex placement by coordinate specification
          \item Now outputs, for convenience, the vertex coordinates of the plotted graph
      }
      \item \code{ridgePathS} has been sped up by killing redundant inversion
      \item The \code{covML} function has been amended with an argument that indicates if a correlation matrix
      (instead of an ML estimate of a covariance matrix) is desired. This offers more flexibility. One can now get the
      ML estimate of the covariance matrix, the ML estimate of the covariance matrix on standardized data, as well as
      the correlation matrix
      \item The \code{optPenalty.LOOCVauto} function has been amended with an argument that indicates if the
      evaluation of the LOOCV score should be performed on the correlation scale
      \item The \code{optPenalty.LOOCV} function has been amended with an argument that indicates if the
      evaluation of the LOOCV score should be performed on the correlation scale
      \item The \code{optPenalty.aLOOCV} function has been amended with an argument that indicates if the
      evaluation of the approximate LOOCV score should be performed on the correlation scale
    }
  }
}


\section{Changes in version 2.0}{
  \subsection{Documentation:}{
    \itemize{
      \item Added this \code{NEWS} file!
      \item Updated (and corrected) \code{CITATION} file
      \item Added \code{README} file
      \item Added (selective) import statements for default packages as required for R-devel
    }
  }
  \subsection{Additions:}{
    \itemize{
      \item \cpkg{rags2ridges} now uses \cpkg{Rcpp} and \cpkg{RcppArmadillo} with
        core functions written in \code{C++}. The package should now be at least
        two orders of magnitude faster in most cases.
      \item Added, next to the core module, the fused ridge module.
        The fused module provides functionality
        for the estimation and graphical modeling of multiple precision matrices from
        multiple high-dimensional data classes.
        Functions from this module are generally suffixed with \code{.fused}.
        Functions tied to (or added with) this module are:
        \itemize{
          \item \code{isSymmetricPD}
          \item \code{isSymmetricPSD}
          \item \code{is.Xlist}
          \item \code{default.target.fused}
          \item \code{createS}
          \item \code{getKEGGPathway}
          \item \code{kegg.target}
          \item \code{pooledS}
          \item \code{pooledP}
          \item \code{KLdiv.fused}
          \item \code{ridgeP.fused}
          \item \code{optPenalty.fused.grid}
          \item \code{print.optPenaltyFusedGrid}
          \item \code{plot.optPenaltyFusedGrid}
          \item \code{optPenalty.fused.auto}
          \item \code{optPenalty.fused}
          \item \code{default.penalty}
          \item \code{fused.test}
          \item \code{print.ptest}
          \item \code{summary.ptest}
          \item \code{hist.ptest}
          \item \code{plot.ptest}
          \item \code{sparsify.fused}
          \item \code{GGMnetworkStats.fused}
          \item \code{GGMpathStats.fused}
        }
      \item The following functions were added to the core module:
        \itemize{
          \item \code{covMLknown}
          \item \code{GGMmutualInfo}
        }
      \item Added miscellaneous (hidden) functions.
    }
  }
  \subsection{Bug fixes:}{
    \itemize{
      \item Fixed bugs in \code{GGMpathstats}:
        \itemize{
          \item Code no longer breaks down if variable names are absent.
          \item Now properly handles singleton pathsets.
        }
      \item Fixed bug in \code{sparsify}: Now always returns symmetric objects
    }
  }
  \subsection{Adjustments:}{
    \itemize{
      \item Argument \code{verticle} as used in various functions has been renamed to
        \code{vertical}. Sorry for any inconvenience.
      \item Internal usage of \code{ridgeS} replaced by the faster C++-dependent counterpart
        \code{ridgeP}
      \item New features updated \code{conditionNumberPlot} function:
        \itemize{
          \item Function has been sped up
          \item Now uses log-equidistant grid for visualization
          \item Now includes option to additionally plot the approximate loss in
            digits of accuracy
        }
    }
  }
  \subsection{Notifications:}{
    \itemize{
      \item \code{ridgeS} is deprecated. Please use \code{ridgeP} instead
      \item Future versions of rags2ridges will be subject to changes in naming conventions
    }
  }
}


\section{Changes in version 1.4}{
  \subsection{Additions:}{
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
    }
  }
  \subsection{Bug fixes:}{
    \itemize{
      \item Corrected small error in \code{evaluateSfit} function.
        The \code{dir} argument was not properly used previously.
    }
  }
  \subsection{Adjustments:}{
    \itemize{
      \item New features updated \code{optPenalty.aLOOCV} function:
        \itemize{
          \item For scalar matrix targets the complete solution path depends on only
            1 eigendecomposition and 1 matrix inversion.
            Meaning: the function is sped up somewhat by lifting redundant
            inversions out of \code{for} loops.
          \item Optional graph now plots the approximated LOOCV negative log-likelihood
            instead of ln(approximated LOOCV negative log-likelihood).
          \item Legend in optional graph has been adapated accordingly.
        }
      \item New features updated \code{optPenalty.LOOCV} function:
        \itemize{
          \item Optional graph now plots the LOOCV negative log-likelihood instead of
            ln(LOOCV negative log-likelihood).
          \item Legend in optional graph has been adapated accordingly.
        }
      \item New features updated \code{default.target} function:
        \itemize{
          \item Inclusion new default target option: \code{type = DIAES}. Gives diagonal
            matrix with inverse of average of eigenvalues of S as entries.
        }
      \item New features updated \code{GGMnetworkStats} function:
        \itemize{
          \item Now also assesses (and returns a logical) if graph/network is chordal.
          \item Now also includes assesment of the eigenvalue centrality.
          \item Now includes option to have list or table output.
        }
      \item New features updated \code{ridgePathS} function:
        \itemize{
          \item Sped up considerably for rotation equivariant alternative estimator.
            By avoidance of redundant eigendecompositions and inversions.
          \item Now catches breakdown due to rounding preculiarities when
            \code{plotType = "pcor"}.
        }
      \item New features updated \code{sparsify} function:
        \itemize{
          \item Inclusion new thresholding function \code{top}: retainment of top elements
            based on absolute partial correlation.
          \item Inclusion output option: When \code{output = "light"}, only the (matrix)
            positions of the zero and non-zero elements are returned.
          \item Function no longer dependent on GeneNet; now makes direct use of
            \cpkg{fdrtool}.
          \item Function now also prints some general information on the number of edges
            retained.
        }
    }
  }
}


\section{Changes in version 1.3}{
  \subsection{Additions:}{
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
    }
  }
  \subsection{Adjustments and name changes:}{
    \itemize{
      \item New features updated \code{evaluateS} function:
        \itemize{
          \item The printed output of the \code{evaluateS} function is now aligned
          \item Calculation spectral condition number has been improved
        }
      \item \code{conditionNumber} function now called \code{conditionNumberPlot}.
        The updated function has new features:
        \itemize{
          \item Main plot can now be obtained with either the spectral (l2) or the
            (approximation to) l1 condition number
          \item Main plot can now be amended with plot of the relative distance to the
            set of singular matrices
          \item The title of the main plot can now be suppressed
          \item One can now obtain numeric output from the function: lambdas and
            condition numbers
        }
      \item New features updated \code{sparsify} function:
        \itemize{
          \item Changed \code{type = c("threshold", "localFDR")} to
            \code{threshold = c("absValue", "localFDR")} (clarifying nomenclature)
          \item Changed \code{threshold} to \code{absValueCut} (clarifying nomenclature)
          \item Will now output both sparsified partial correlation/standardized
            precision matrix and the sparsified precison matrix,
            when input consists of the unstandardized precision matrix
        }
      \item New features updated \code{ridgeS} function:
        \itemize{
          \item Contains an improved evaluation of the target matrix possibly being a
            null matrix
          \item Now evaluates if a rotation equivariant alternative estimator ensues for
            a given target matrix
          \item When rotation equivariant alternative estimator ensues, computation is
            sped up considerably by circumventing the matrix square root
        }
      \item \code{optPenaltyCV} function now called \code{optPenalty.LOOCV}, for sake
        of (naming) consistency. The updated function has new features:
        \itemize{
          \item \code{targetScale} option has been removed
          \item Replaced \code{log} in optional graph by \code{ln}
          \item Visual layout of optional graph now more in line with recommendations by
            Tufte (regarding data-ink ratio)
        }
      \item New features updated \code{optPenalty.aLOOCV} function:
        \itemize{
          \item Replaced \code{log} in optional graph by \code{ln}
          \item Visual layout of optional graph now more in line with recommendations by
            Tufte (regarding data-ink ratio)
        }
      \item Computation optimal penalty in \code{conditionNumberPlot},
        \code{optPenalty.aLOOCV} and \code{optPenalty.LOOCV} functions sped up
        considerably for rotation equivariant alternative estimator.
        By usage new ridgeS and avoidance of redundant eigendecompositions
      \item Default target in \code{ridgeS}, \code{conditionNumberPlot},
        \code{optPenalty.aLOOCV} and \code{optPenalty.LOOCV} now \code{"DAIE"}
        option from \code{default.target}
    }
  }
}


\section{Changes in version 1.2}{
  \subsection{Additions:}{
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
    }
  }
  \subsection{Bug fixes:}{
    \itemize{
      \item Corrected small error in \code{optPenaltyCV} function
    }
  }
  \subsection{Adjustments:}{
    \itemize{
      \item Both \code{optPenaltyCV} and \code{optPenalty.aLOOCV} now utilize
        \code{covML} instead of \code{cov}
      \item Default output option in \code{optPenaltyCV} (as in
        \code{optPenalty.aLOOCV}) is now \code{light}
    }
  }
}
