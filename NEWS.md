# News for Package [rags2ridges](https://cran.r-project.org/package=rags2ridges)


# Version 2.2.7

## Adjustments:
  * Functions optPenaltyPchordal, ridgePchordal, ridgePsign, and support4ridgeP (temporarily) deprecated (for major adjustments)
  * Replaced if() conditions comparing class() to string with evaluations using inherits()
  
## Documentation:
  * Fixed URLs.


# Version 2.2.6

## Documentation:
  * Canonicalization of URLs.
  * Update of published papers 


# Version 2.2.5

## Documentation:
  * Improved documentation and added new [pkgdown](https://cran.r-project.org/package=pkgdown) documentation website.
  * NEWS file moved to markdown format instead of .Rd and available on the website


# Version 2.2.4

## Adjustments:
  *  Documentation roxygenized.
  *  More selective importing and exporting.
  *  S3 implementation of `ridgeP` output.


# Version 2.2.3

## Documentation:
  * Updated `CITATION` file
  * Updated `README` file

## Bug fixes:
  * Fixed bug in `GGMpathStats`:
    Incorrectly stated before that all igraph layouts were supported.
    Now they indeed are supported.

## Adjustments:
  * Bioconductor dependencies are now automatically installed upon first 
    installation of rags2ridges.
  * `GGMpathStats` now has additional visualization options: It can handle all 
    layout functions supported by igraph. Moreover, it is now possible to 
    specify custom coordinates for node-placement.


# Version 2.2.2

## Notifications:
  * Hot fix due to class changes in "matrix". No major visible user changes.
  * `CNplot` function again updated: higher max. iterations for Lanczos method

# Version 2.2.1

## Notifications:
  * Hot fix due to new RNG. No visible user changes.


# Version 2.2

## Notifications:
  * `optPenalty.LOOCV` is deprecated. Please use `optPenalty.kCV` instead
  * `optPenalty.LOOCVauto` is deprecated. Please use `optPenalty.kCVauto` instead


# Version 2.1.1

## Documentation:
  * Updated `CITATION` file
  * Updated `README` file

## Adjustments:
  * `sparsify` now has an additional thresholding option: 'connected'


# Version 2.1

## Documentation:
  * Updated `CITATION` file
  * Updated `README` file

## Bug fixes:
  * Fixed bug in `Ugraph`:
    - Incorrectly stated before that all igraph layouts were supported.
    - Now they indeed are supported.

## Notifications:
  * `conditionNumberPlot` is deprecated. Please use `CNplot` instead
  * Features of the `CNplot` function (above and beyond `conditionNumberPlot`):
    - The `digitLoss` and `rlDist` arguments have been removed
    - These arguments have been replaced with the logical argument `Iaids`
    - `Iaids = TRUE` amends the basic condition number plot with interpretational aids
    - These aids are the approximate loss in digits of accuracy and and approximation of the acceleration             along the regularization path of the condition number
    - Argument `main` is now a character argument
    - Argument `value` now by default takes the value 1e-100 (convenient)
    - Now uses C++ functionalty for additional speed

## Adjustments:
  * `edgeHeat` now has aligned x-axis labels
  * The visualizations of the `optPenalty.LOOCV` and `optPenalty.aLOOCV` functions
      now will no longer produce horizontal and/or vertical lines that fall outside the boundaries
      of the figure
  * `optPenalty.LOOCV` now uses log-equidistant penalty grid for optimal penalty
      parameter determination (this also enhances the visualization)
  * New features updated `optPenalty.aLOOCV` function:
    - Function has been sped up by killing redundant inversion
    - now uses log-equidistant penalty grid for optimal penalty
          parameter determination (this also enhances the visualization)
  * New features updated `Ugraph` function:
    - One can now also specify vertex placement by coordinate specification
    - Now outputs, for convenience, the vertex coordinates of the plotted graph
  * `ridgePathS` has been sped up by killing redundant inversion
  * The `covML` function has been amended with an argument that indicates if a correlation matrix
      (instead of an ML estimate of a covariance matrix) is desired. This offers more flexibility. One can now get the
      ML estimate of the covariance matrix, the ML estimate of the covariance matrix on standardized data, as well as
      the correlation matrix
  * The `optPenalty.LOOCVauto` function has been amended with an argument that indicates if the
      evaluation of the LOOCV score should be performed on the correlation scale
  * The `optPenalty.LOOCV` function has been amended with an argument that indicates if the
      evaluation of the LOOCV score should be performed on the correlation scale
  * The `optPenalty.aLOOCV` function has been amended with an argument that indicates if the
      evaluation of the approximate LOOCV score should be performed on the correlation scale


# Version 2.0

## Documentation:
  * Added this `NEWS` file!
  * Updated (and corrected) `CITATION` file
  * Added `README` file
  * Added (selective) import statements for default packages as required for R-devel

## Additions:
  * [rags2ridges](https://cran.r-project.org/package=rags2ridges) now uses [Rcpp](https://cran.r-project.org/package=Rcpp) and [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo) with
    core functions written in `C++`. The package should now be at least
    two orders of magnitude faster in most cases.
  * Added, next to the core module, the fused ridge module.
    The fused module provides functionality
    for the estimation and graphical modeling of multiple precision matrices from
    multiple high-dimensional data classes.
    Functions from this module are generally suffixed with `.fused`.
    Functions tied to (or added with) this module are:
    - `isSymmetricPD`
    - `isSymmetricPSD`
    - `is.Xlist`
    - `default.target.fused`
    - `createS`
    - `getKEGGPathway`
    - `kegg.target`
    - `pooledS`
    - `pooledP`
    - `KLdiv.fused`
    - `ridgeP.fused`
    - `optPenalty.fused.grid`
    - `print.optPenaltyFusedGrid`
    - `plot.optPenaltyFusedGrid`
    - `optPenalty.fused.auto`
    - `optPenalty.fused`
    - `default.penalty`
    - `fused.test`
    - `print.ptest`
    - `summary.ptest`
    - `hist.ptest`
    - `plot.ptest`
    - `sparsify.fused`
    - `GGMnetworkStats.fused`
    - `GGMpathStats.fused`
  * The following functions were added to the core module:
    - `covMLknown`
    - `GGMmutualInfo`
  * Added miscellaneous (hidden) functions.

## Bug fixes:
  * Fixed bugs in `GGMpathstats`:
    - Code no longer breaks down if variable names are absent.
    - Now properly handles singleton pathsets.
  * Fixed bug in `sparsify`: Now always returns symmetric objects

## Adjustments:
  * Argument `verticle` as used in various functions has been renamed to
    `vertical`. Sorry for any inconvenience.
  * Internal usage of `ridgeS` replaced by the faster C++-dependent counterpart
    `ridgeP`
  * New features updated `conditionNumberPlot` function:
    - Function has been sped up
    - Now uses log-equidistant grid for visualization
    - Now includes option to additionally plot the approximate loss in
      digits of accuracy

## Notifications:
  * `ridgeS` is deprecated. Please use `ridgeP` instead
  * Future versions of rags2ridges will be subject to changes in naming conventions


# Version 1.4

## Additions:
  * Inclusion hidden function `.pathContribution` for usage in
    `GGMpathStats` function
  * Inclusion hidden function `.path2string` for usage in
    `GGMpathStats` function
  * Inclusion hidden function `.pathAndStats` for usage in
    `GGMpathStats` function
  * Inclusion hidden function `.cvl` for usage in
    `optPenalty.LOOCVauto` function
  * Inclusion hidden function `.lambdaNullDist` for usage in
    `GGMblockNullPenalty` function
  * Inclusion hidden function `.blockTestStat` for usage in
    `GGMblockTest` function
  * Inclusion function that expresses the covariance between a pair of
    variables as a sum of path weights: `GGMpathStats`
  * Inclusion function that determines the optimal penalty parameter
    value by application of the Brent algorithm to the LOOCV log-likelihood:
    `optPenalty.LOOCVauto`
  * Inclusion function that generates the distribution of the penalty
    parameter under the null hypothesis of block independence:
    `GGMblockNullPenalty`
  * Inclusion function that performs a permutation test for block
    structure in the precision matrix: `GGMblockTest`
  * Inclusion wrapper function: `fullMontyS`

## Bug fixes:
  * Corrected small error in `evaluateSfit` function.
    The `dir` argument was not properly used previously.

## Adjustments:
  * New features updated `optPenalty.aLOOCV` function:
    - For scalar matrix targets the complete solution path depends on only
      1 eigendecomposition and 1 matrix inversion.
      Meaning: the function is sped up somewhat by lifting redundant
      inversions out of `for` loops.
    - Optional graph now plots the approximated LOOCV negative log-likelihood
      instead of ln(approximated LOOCV negative log-likelihood).
    - Legend in optional graph has been adapated accordingly.
  * New features updated `optPenalty.LOOCV` function:
    - Optional graph now plots the LOOCV negative log-likelihood instead of
      ln(LOOCV negative log-likelihood).
    - Legend in optional graph has been adapated accordingly.
  * New features updated `default.target` function:
    - Inclusion new default target option: `type = DIAES`. Gives diagonal
      matrix with inverse of average of eigenvalues of S as entries.
  * New features updated `GGMnetworkStats` function:
    - Now also assesses (and returns a logical) if graph/network is chordal.
    - Now also includes assesment of the eigenvalue centrality.
    - Now includes option to have list or table output.
  * New features updated `ridgePathS` function:
    - Sped up considerably for rotation equivariant alternative estimator.
      By avoidance of redundant eigendecompositions and inversions.
    - Now catches breakdown due to rounding preculiarities when
      `plotType = "pcor"`.
  * New features updated `sparsify` function:
    - Inclusion new thresholding function `top`: retainment of top elements
      based on absolute partial correlation.
    - Inclusion output option: When `output = "light"`, only the (matrix)
      positions of the zero and non-zero elements are returned.
    - Function no longer dependent on GeneNet; now makes direct use of
      [fdrtool](https://cran.r-project.org/package=fdrtool).
    - Function now also prints some general information on the number of edges
      retained.


# Version 1.3

## Additions:
  * Inclusion hidden function `.ridgeSi` for usage in
    `conditionNumberPlot` function.
  * Inclusion hidden function `.eigShrink` for usage in (a.o.)
    `ridgeS` function.
  * Inclusion function calculating various network statistics from a
    sparse matrix: `GGMnetworkStats`
  * Inclusion function for visual inspection fit of regularized precision
    matrix to sample covariance matrix: `evaluateSfit`
  * Inclusion function for visualization of regularization paths:
    `ridgePathS`
  * Inclusion function for default target matrix generation:
    `default.target`

## Adjustments and name changes:
  * New features updated `evaluateS` function:
    - The printed output of the `evaluateS` function is now aligned
    - Calculation spectral condition number has been improved
  * `conditionNumber` function now called `conditionNumberPlot`.
    The updated function has new features:
    - Main plot can now be obtained with either the spectral (l2) or the
      (approximation to) l1 condition number
    - Main plot can now be amended with plot of the relative distance to the
      set of singular matrices
    - The title of the main plot can now be suppressed
    - One can now obtain numeric output from the function: lambdas and
      condition numbers
  * New features updated `sparsify` function:
    - Changed `type = c("threshold", "localFDR")` to
      `threshold = c("absValue", "localFDR")` (clarifying nomenclature)
    - Changed `threshold` to `absValueCut` (clarifying nomenclature)
    - Will now output both sparsified partial correlation/standardized
      precision matrix and the sparsified precison matrix,
      when input consists of the unstandardized precision matrix
  * New features updated `ridgeS` function:
    - Contains an improved evaluation of the target matrix possibly being a
      null matrix
    - Now evaluates if a rotation equivariant alternative estimator ensues for
      a given target matrix
    - When rotation equivariant alternative estimator ensues, computation is
      sped up considerably by circumventing the matrix square root
  * `optPenaltyCV` function now called `optPenalty.LOOCV`, for sake
    of (naming) consistency. The updated function has new features:
    - `targetScale` option has been removed
    - Replaced `log` in optional graph by `ln`
    - Visual layout of optional graph now more in line with recommendations by
      Tufte (regarding data-ink ratio)
  * New features updated `optPenalty.aLOOCV` function:
    - Replaced `log` in optional graph by `ln`
    - Visual layout of optional graph now more in line with recommendations by
      Tufte (regarding data-ink ratio)
  * Computation optimal penalty in `conditionNumberPlot`,
    `optPenalty.aLOOCV` and `optPenalty.LOOCV` functions sped up
    considerably for rotation equivariant alternative estimator.
    By usage new ridgeS and avoidance of redundant eigendecompositions
  * Default target in `ridgeS`, `conditionNumberPlot`,
    `optPenalty.aLOOCV` and `optPenalty.LOOCV` now \code{"DAIE"
    option from `default.target`


# Version 1.2

## Additions:
  * Inclusion function for ML estimation of the sample covariance matrix:
    `covML`
  * Inclusion function for approximate leave-one-out cross-validation:
    `optPenalty.aLOOCV`
  * Inclusion function `conditionNumber` to visualize the
    spectral condition number over the regularization path
  * Inclusion function `evaluateS` to evaluate basic
    properties of a covariance matrix
  * Inclusion function `KLdiv` that calculates the
    Kullback-Leibler divergence between two normal distributions
  * Inclusion option to suppress on-screen output in `sparsify`
    function

## Bug fixes:
  * Corrected small error in `optPenaltyCV` function

## Adjustments:
  * Both `optPenaltyCV` and `optPenalty.aLOOCV` now utilize
    `covML` instead of `cov`
  * Default output option in `optPenaltyCV` (as in
    `optPenalty.aLOOCV`) is now `light`
