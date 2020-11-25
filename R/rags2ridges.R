################################################################################
################################################################################
################################################################################
##
## Name: rags2ridges
## Authors: Carel F.W. Peeters, Anders E. Bilgrau, and Wessel N. van Wieringen
##
## Maintainer: Carel F.W. Peeters
##             Statistics for Omics Research Unit
##             Dept. of Epidemiology & Biostatistics
##             Amsterdam Public Health research institute
##             VU University medical center
##             Amsterdam, the Netherlands
## Email:	     cf.peeters@vumc.nl
##
## Version: 2.2.3
## Last Update:	18/08/2020
## Description:	Ridge estimation for high-dimensional precision matrices
##              Includes supporting functions for (integrative) graphical modeling
##
## Code files: rags2ridges.R      >> master file/core module
##             rags2ridgesFused.R >> fused module
##             rags2ridgesMisc.R  >> miscellaneous module
##             rags2ridgesDepr.R  >> deprecated function collection
##             rags2ridges.cpp    >> C++ work horses
##
## Publications:
##   [1] van Wieringen, W.N. & Peeters, C.F.W. (2016).
##       "Ridge Estimation of Inverse Covariance Matrices from High-Dimensional
##       Data", Computational Statistics & Data Analysis, vol 103: 284-303.
##   [2] van Wieringen, W.N. & Peeters, C.F.W. (2015).
##       "Application of a New Ridge Estimator of the Inverse Covariance Matrix
##       to the Reconstruction of Gene-Gene Interaction Networks", in: di Serio,
##       C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.) Computational
##       Intelligence Methods for Bioinformatics and Biostatistics. Lecture
##       Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
##   [3] Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., &
##       van Wieringen, W.N. (2020).
##       "Targeted Fused Ridge Estimation of Inverse Covariance Matrices from
##       Multiple High-Dimensional Data Classes", Journal of Machine Learning
##       Research, vol. 21(26): 1-52.
## 	 [4] Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2020).
##       "The Spectral Condition Number Plot for Regularization Parameter
##       Evaluation", Computational Statistics, vol. 35: 629-646.
##
################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module A: rags2ridges Core
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

##------------------------------------------------------------------------------
##
## Hidden support functions
##
##------------------------------------------------------------------------------

.trace <- function(M){
  ##############################################################################
  # - Internal function to compute the trace of a matrix
  # - Faster support function (as opposed to 'matrix.trace') when input M is
  #   already forced to 'matrix'
  # - M > matrix input
  ##############################################################################

  return(sum(diag(M)))
}



.is.int <- function(x, tolerance = .Machine$double.eps){
  ##############################################################################
  # - Logical function that checks if a number is an integer within machine
  #   precision
  # - x         > input number
  # - tolerance > tolerance threshold for determining integer quality
  ##############################################################################

  #return(all.equal(x, as.integer(x), ...)) # Faster(?), safer
  return(abs(x - round(x)) < tolerance)
}



.LL <- function(S, P){
  ##############################################################################
  # - Function that computes the value of the (negative) log-likelihood
  # - S > sample covariance matrix
  # - P > precision matrix (possibly regularized inverse of covariance or
  #       correlation matrix)
  ##############################################################################

  LL <- -log(det(P)) + sum(S*P) #.trace(S %*% P)
  return(LL)
}



.Frobenius <- function(X) {
  ##############################################################################
  # - Function computing squared Frobenius norm - the sum of the squarred
  #   entries.
  # - X > A numeric
  ##############################################################################

  return(sum(X^2))
}



.FrobeniusLoss <- function(O, P){
  ##############################################################################
  # - Function computing Frobenius loss
  # - O > Estimated (possibly regularized) precision matrix
  # - P > True (population) precision matrix
  ##############################################################################

  return(.Frobenius(O - P))
}



.RelativeFrobeniusLoss <- function(O, P) {
  ##############################################################################
  # - Function computing the frobenius loss relative to the (frobenious) norm
  #   of P.
  # - O > Estimated (possibly regularized) precision matrix
  # - P > True (population) precision matrix
  ##############################################################################

  return(.FrobeniusLoss(O, P)/.Frobenius(P))
}



.QuadraticLoss <- function(O, C){
  ##############################################################################
  # - Function computing Quadratic loss
  # - O > Estimated (possibly regularized) precision matrix
  # - C > True (population) covariance matrix
  ##############################################################################

  return((sum(((O %*% C - diag(ncol(O))))^2)))
}



.eigShrink <- function(dVec, lambda, const = 0){
  ##############################################################################
  # - Function that shrinks the eigenvalues in an eigenvector
  # - Shrinkage is that of rotation equivariant alternative ridge estimator
  # - Main use is in avoiding expensive matrix square root when choosing a
  #   target that leads to a rotation equivariant version of the alternative
  #   ridge estimator
  # - dVec   > numeric vector containing the eigenvalues of a matrix S
  # - lambda > penalty parameter
  # - const  > a constant, default = 0
  ##############################################################################

  Evector <- (dVec - lambda * const)
  return(sqrt(lambda + Evector^2/4) + Evector/2)
}



.ridgeSi <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Hidden function that calculates Ridge estimators of a covariance matrix
  # - Function is mirror image main routine 'ridgeS'
  # - Main use is to circumvent (unnecessary) inversion (especially in
  #   'conditionNumberPlot' function)
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a
  #             covariance matrix
  # - ArchI   > Archetypal I ridge estimator of a covariance matrix
  # - ArchII  > Archetypal II ridge estimator of a covariance matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II estimator
  # - When target is not the null-matrix it is expected to be p.d. for the vWP
  #   type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  # Alternative estimator
  if (type == "Alt"){
    if (all(target == 0)){
      Spectral  <- eigen(S, symmetric = TRUE)
      Eigshrink <- .eigShrink(Spectral$values, lambda)
      C_Alt     <- Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)
      colnames(C_Alt) = rownames(C_Alt) <- colnames(S)
      return(C_Alt)
    } else if (all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
      varPhi    <- unique(diag(target))
      Spectral  <- eigen(S, symmetric = TRUE)
      Eigshrink <- .eigShrink(Spectral$values, lambda, const = varPhi)
      C_Alt     <- Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)
      colnames(C_Alt) = rownames(C_Alt) <- colnames(S)
      return(C_Alt)
    } else {
      D     <- (S - lambda * target)
      C_Alt <- D/2 + sqrtm((D %*% D)/4 + lambda * diag(nrow(S)))
      return(C_Alt)
    }
  }

  # Archetypal I
  if (type == "ArchI"){
    C_ArchI <- (1-lambda) * S + lambda * solve(target)
    return(C_ArchI)
  }

  # Archetypal II
  if (type == "ArchII"){
    C_ArchII <- S + lambda * diag(nrow(S))
    return(C_ArchII)
  }
}



.pathContribution <- function(sparseP, path, detSparseP){
  ##############################################################################
  # - Function calculating the contribution of a path to the covariance between
  #   begin and end node
  # - sparseP    > sparse precision/partial correlation matrix
  # - path       > path between two nodes (start and end node)
  # - detSparseP > determinant of 'sparseP'
  ##############################################################################

  if (length(path) < nrow(sparseP)){
    return((-1)^(1+length(path)) *
             prod(sparseP[cbind(path[-1], path[-length(path)])]) *
             det(sparseP[-path, -path]) / detSparseP)
  }
  if (length(path) == nrow(sparseP)){
    return((-1)^(1+length(path)) *
             prod(sparseP[cbind(path[-1], path[-length(path)])]) / detSparseP)
  }
}



.path2string <- function(path){
  ##############################################################################
  # - Function converting a numeric or character vector into a single string
  # - path > path between two nodes (start and end node)
  ##############################################################################

  pName <- sprintf("%s--%s", path[1], path[2])
  if (length(path) > 2){
    for (w in 3:(length(path))){ pName <- sprintf("%s--%s", pName, path[w]) }
  }
  return(pName)
}



.pathAndStats <- function(Gt, node1t, node2t, nei1t, nei2t, P0t, detP0t,
                          pathNames){
  ##############################################################################
  # - Function determining shortest paths (and their contribution) between
  #   node 1 and node 2
  # - It does so via the neighborhoods 1 and 2
  # - Gt        > graphical object
  # - node1t    > start node of the path
  # - node2t    > end node of the path
  # - nei1t     > neighborhood around the start node
  # - nei2t     > neighborhood around the end node
  # - p0t       > sparse precision/partial correlation matrix
  # - detP0t    > determinant of 'p0t'
  # - pathNames > named path represented as string
  ##############################################################################

  pathsTemp <- list()
  pathStatsTemp <- numeric()
  for (v1 in 1:length(nei1t)){
    for (v2 in 1:length(nei2t)){
      slhNo1Ne1 <- get.all.shortest.paths(Gt, node1t, nei1t[v1])$res
      slhNe1Ne2 <- get.all.shortest.paths(Gt, nei1t[v1], nei2t[v2])$res
      slhNe2No2 <- get.all.shortest.paths(Gt, nei2t[v2], node2t)$res
      for (uNo1Ne1 in 1:length(slhNo1Ne1)){
        for (uNe1Ne2 in 1:length(slhNe1Ne2)){
          for (uNe2No2 in 1:length(slhNe2No2)){
            fullPath <- c(slhNo1Ne1[[uNo1Ne1]], slhNe1Ne2[[uNe1Ne2]][-1],
                          slhNe2No2[[uNe2No2]][-1])
            if (length(unique(fullPath)) == (length(fullPath))){
              pName <- .path2string(fullPath)
              if (!(pName %in% c(pathNames, rownames(pathStatsTemp)))){
                pathsTemp[[length(pathsTemp)+1]] <- fullPath
                pathStatsTemp <-
                  rbind(pathStatsTemp,
                        c(length(fullPath)-1,
                          .pathContribution(P0t, fullPath, detP0t)))
                rownames(pathStatsTemp)[nrow(pathStatsTemp)] <- pName
              }
            }
          }
        }
      }
    }
  }
  return(list(paths=pathsTemp, pathStats=pathStatsTemp))
}



.cvl <- function(lambda, Y, cor = FALSE, target = default.target(covML(Y)),
                 type = "Alt"){
  ##############################################################################
  # - Function that calculates a cross-validated negative log-likelihood score
  #   for single penalty value
  # - lambda > value penalty parameter
  # - Y      > (raw) Data matrix, variables in columns
  # - cor    > logical indicating if evaluation of the LOOCV score should be
  #            performed on the correlation matrix
  # - target > target (precision terms) for Type I estimators,
  #            default = default.target(covML(Y))
  # - type   > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  ##############################################################################

  slh <- numeric()
  for (i in 1:nrow(Y)){
    S   <- covML(Y[-i, ], cor = cor)
    slh <- c(slh, .LL(t(Y[i, , drop = FALSE]) %*% Y[i, , drop = FALSE],
                      ridgeP(S, lambda, target = target, type = type)))
  }
  return(mean(slh))
}



.kcvl <- function(lambda, Y, cor, target, type, folds){
  ##############################################################################
  # - Function that calculates a cross-validated negative log-likelihood score
  #   for single penalty value
  # - lambda > value penalty parameter
  # - Y      > (raw) Data matrix, variables in columns
  # - cor    > logical indicating if evaluation of the LOOCV score should be
  #            performed on the correlation matrix
  # - target > target (precision terms) for Type I estimators,
  #            default = default.target(covML(Y))
  # - type   > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - folds  > cross-validation sample splits
  ##############################################################################

  cvLL <- 0
  for (f in 1:length(folds)){
    S   <- covML(Y[-folds[[f]], , drop=FALSE], cor = cor)
    cvLL <- cvLL + .LL(crossprod(Y[folds[[f]], , drop=FALSE]) / length(folds[[f]]),
                       ridgeP(S, lambda, target=target, type=type))
  }
  return(cvLL / length(folds))
}



.lambdaNullDist <- function(i, Y, id, lambdaMin, lambdaMax,
                            lambdaInit, target, type){
  ##############################################################################
  # - Function that determines the optimal value of the penalty parameter for a
  #   single permutation
  # - Optimal penalty determined using the 'optPenalty.LOOCVauto' function
  # - i          > number of permutations; passed by nPerm in
  #                'GGMblockNullPenalty'
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - target     > target (precision terms) for Type I estimators
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}
  ##############################################################################

  reshuffle    <- sample(1:nrow(Y), nrow(Y))
  Y[, id == 1] <- Y[reshuffle, id == 1]
  return(optPenalty.LOOCVauto(Y, lambdaMin, lambdaMax, lambdaInit,
                              target = target, type = type)$optLambda)
}



.blockTestStat <- function(i, Y, id, lambda, target, type){
  ##############################################################################
  # - Function that calculates logratio statistic for block independence
  # - i      > number of permutations; passed by nPerm in 'GGMblockTest'
  # - Y      > (raw) Data matrix, variables in columns
  # - id     > indicator variable for the two blocks of the precision matrix
  # - lambda > value penalty parameter
  # - target > target (precision terms) for Type I estimators
  # - type   > must be one of {"Alt", "ArchI", "ArchII"}
  ##############################################################################

  reshuffle    <- sample(1:nrow(Y), nrow(Y))
  Y[, id == 1] <- Y[reshuffle, id == 1]
  S <- solve(ridgeP(covML(Y), lambda = lambda, target = target, type = type))
  return(log(det(S[id == 0, id == 0])) +
           log(det(S[id == 1, id == 1])) - log(det(S)))
}




##------------------------------------------------------------------------------
##
## Support functions
##
##------------------------------------------------------------------------------







#' Symmetrize matrix
#' 
#' Function that symmetrizes matrices.
#' 
#' Large objects that are symmetric sometimes fail to be recognized as such by
#' R due to rounding under machine precision. This function symmetrizes for
#' computational purposes matrices that are symmetric in numeric ideality.
#' 
#' @param M (In numeric ideality symmetric) square \code{matrix}.
#' @return A symmetric \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCV(X, 10, 30, 10, target = diag(diag(1/covML(X))))
#' 
#' ## Check symmetry
#' ## OPT$optPrec is symmetric by definition
#' ## But is not recognized as such due to rounding peculiarities
#' isSymmetric(OPT$optPrec)
#' 
#' ## Symmetrize
#' symm(OPT$optPrec)
#' 
#' @export symm
symm <- function(M){
  ##############################################################################
  # - Large objects that are symmetric sometimes fail to be recognized as such
  #   by R due to rounding under machine precision. This function symmetrizes
  #   for computational purposes matrices that are symmetric in numeric ideality
  # - M > symmetric (in numeric ideality) square matrix
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be a square matrix")
  }
  else {
    # Symmetrize
    Msym <- (M + t(M))/2

    # Return
    return(Msym)
  }
}









#' Transform real matrix into an adjacency matrix
#' 
#' Function that transforms a real matrix into an adjacency matrix. Intended
#' use: Turn sparsified precision matrix into an adjacency matrix for
#' undirected graphical representation.
#' 
#' 
#' @param M (Possibly sparsified precision) \code{matrix}.
#' @param diag A \code{logical} indicating if the diagonal elements should be
#' retained.
#' @return Function returns an adjacency \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}}, \code{\link{sparsify}},
#' \code{\link{edgeHeat}}, \code{\link{Ugraph}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain regularized precision matrix
#' P <- ridgeP(Cx, lambda = 10, type = "Alt")
#' 
#' ## Obtain sparsified partial correlation matrix
#' PC0 <- sparsify(P, threshold = "localFDR", FDRcut = .8)
#' 
#' ## Obtain adjacency matrix
#' adjacentMat(PC0$sparsePrecision)
#' 
#' @export adjacentMat
adjacentMat <- function(M, diag = FALSE){
  ##############################################################################
  # - Function that transforms a real matrix into an adjacency matrix
  # - Intended use: Turn sparsified precision matrix into an adjacency matrix
  #   for undirected graph
  # - M    > (sparsified precision) matrix
  # - diag > logical indicating if the diagonal elements should be retained
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be square matrix")
  }
  else {
    # Create adjacency matrix
    AM <- M
    AM[AM != 0] <- 1
    diag(AM) <- 0

    if (diag){
      diag(AM) <- 1
    }

    # Return
    return(AM)
  }
}









#' Maximum likelihood estimation of the covariance matrix
#' 
#' Function that gives the maximum likelihood estimate of the covariance
#' matrix.
#' 
#' The function gives the maximum likelihood (ML) estimate of the covariance
#' matrix. The input matrix \code{Y} assumes that the variables are represented
#' by the columns. Note that when the input data is standardized, the ML
#' covariance matrix of the scaled data is computed. If a correlation matrix is
#' desired, use \code{cor = TRUE}.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param cor A \code{logical} indicating if the correlation matrix should be
#' returned
#' @return Function returns the maximum likelihood estimate of the covariance
#' \code{matrix}. In case \code{cor = TRUE}, the correlation matrix is
#' returned.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain ML estimate covariance matrix
#' Cx <- covML(X)
#' 
#' ## Obtain correlation matrix
#' Cx <- covML(X, cor = TRUE)
#' 
#' @export covML
covML <- function(Y, cor = FALSE){
  ##############################################################################
  # - function that gives the maximum likelihood estimate of the covariance
  # - Y   > (raw) data matrix, assumed to have variables in columns
  # - cor > logical indicating if the correlation matrix should be returned
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (class(cor) != "logical"){
    stop("Input (cor) is of wrong class")
  }
  else {
    if (cor){
      Sml <- cor(Y)
    }
    else {
      Ys  <- scale(Y, center = TRUE, scale = FALSE)
      Sml <- crossprod(Ys)/nrow(Ys)  # (t(Ys) %*% Ys)/nrow(Ys)
    }

    # Return
    return(Sml)
  }
}









#' Maximum likelihood estimation of the covariance matrix with assumptions on
#' its structure
#' 
#' Function that performs maximum likelihood estimation of the covariance
#' matrix, with various types of assumptions on its structure.
#' 
#' The function gives the maximum likelihood estimate of the covariance matrix.
#' The input matrix \code{Y} assumes that the variables are represented by the
#' columns.
#' 
#' When simultaneously \code{covMat=NULL}, \code{corMat=NULL},
#' \code{corType="none"} and \code{varType="none"} the \code{covML}-function is
#' invoked and the regular maximum likelihood estimate of the covariance matrix
#' is returned.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param covMat A positive-definite covariance \code{matrix}. When specified,
#' the to-be-estimated covariance matrix is assumed to be proportional to the
#' specified covariance matrix. Hence, only a constant needs to estimated.
#' @param corMat A positive-definite correlation \code{matrix}. When specified,
#' the to-be-estimated covariance matrix is assumed to have this correlation
#' structure. Hence, only the marginal variances need to be estimated.
#' @param corType A \code{character}, either \code{"none"} (no structure on the
#' correlation among variate assumed) or \code{"equi"} (variates are
#' equi-correlated).
#' @param varType A \code{character}, either \code{"none"} (no structure on the
#' marginal variances of the variates assumed) or \code{"common"} (variates
#' have equal marginal variances).
#' @param nInit An \code{integer} specifying the maximum number of iterations
#' for likelihood maximization when \code{corType="equi"} .
#' @return The maximum likelihood estimate of the covariance \code{matrix}
#' under the specified assumptions on its structure.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{covML}}
#' @examples
#' 
#' ## Obtain some data
#' p = 10
#' n = 100
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:10] = letters[1:10]
#' 
#' ## Obtain maximum likelihood estimate covariance matrix
#' Cx <- covMLknown(X, corType="equi", varType="common")
#' 
#' @export covMLknown
covMLknown <- function(Y, covMat = NULL, corMat = NULL,
                       corType = "none", varType = "none", nInit = 100){
  ##############################################################################
  # - Maximum likelihood estimation of the covariance matrix, with various
  #   types of assumptions on the structure of this matrix.
  # - Y	      > (raw) data matrix, assumed to have variables in columns
  # - covMat  > A positive-definite covariance 'matrix'. When specified, the
  #             to-be-estimated covariance matrix is assumed to be proportional
  #             to the specified covariance matrix. Hence, only a constant needs
  #             to be estimated
  # - corMat  > A positive-definite correlation 'matrix'. When specified, the
  #             to-be-estimated covariance matrix is assumed to have this
  #             correlation structure. Hence, only the variances need to
  #		          be estimated.
  # - corType > character that is either "none" (no structure on the correlation
  #             among variate assumed) or "equi" (variates are equi-correlated).
  # - varType	> character that is either "none" (no structure on the marginal
  #             variances of the variates assumed) or "common" (variates have
  #             equal marginal variances).
  # - nInit   > Maximum number of iterations for likelihood maximization
  #             when corType='equi'.
  #
  # NOTES:
  # - Future version should a.o. also allow a first order autoregressive
  #   correlation assumption.
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  # center data
  Ys <- scale(Y, center = TRUE, scale = FALSE)

  # equicorrelated covariance: parameter estimation
  if (corType=="equi" & varType !="none"){

    # initial variance estimate
    sds <- sqrt(apply(Ys^2, 2, mean))
    rho <- 0
    llNew <- 10^(-10)
    for (k in 1:nInit){
      # update current log-likelihood
      llPrev <- llNew

      # estimate rho
      Sml <- diag(1/sds) %*% ((t(Ys) %*% Ys)/nrow(Ys)) %*% diag(1/sds)
      a1 <- sum(diag(Sml))
      a2 <- (ncol(Ys)-1) * a1 - sum(Sml)
      p <- nrow(Sml)
      minLoglikEqui <- function(rho, p, a1, a2){
        (p-1) * log(1-rho) +  log((p-1) * rho + 1) + ((1-rho) * ((p-1) * rho + 1))^(-1) * (a1 + a2 * rho) }
      rho <- optim(par=0.1, fn=minLoglikEqui, method="Brent", lower=-1/(p-1)
                   + .Machine$double.eps, upper=1-.Machine$double.eps, p=ncol(Ys),
                   a1=a1, a2=a2)$par
      llNew <- minLoglikEqui(rho, p, a1, a2)
      if (abs(llNew - llPrev) < 0.0001){ break }

      # estimate variance(s)
      Sml <- matrix(rho, ncol(Ys), ncol(Ys))
      diag(Sml) <- 1
      if (varType!="common"){
        V <- eigen(Sml)
        D <- abs(V$values)
        V <- V$vectors
        Vinner <- eigen(diag(1/sqrt(D)) %*% t(V) %*% ((t(Ys) %*% Ys)/nrow(Ys))
                        %*% V %*% diag(1/sqrt(D)))
        Vinner$values <- sqrt(abs(Re(Vinner$values)))
        sds <- Re(diag(V %*% diag(sqrt(D)) %*% Vinner$vectors
                       %*% diag(Vinner$values) %*% t(Vinner$vectors)
                       %*% diag(sqrt(D)) %*% t(V)))
      }
      if (varType=="common"){
        sds <- rep(sum(diag(solve(Sml) %*% (t(Ys) %*% Ys)/nrow(Ys))) /
                     ncol(Ys), ncol(Ys))
      }
    }

    # equicorrelated covariance: matrix construction
    Sml <- matrix(rho, ncol(Ys), ncol(Ys))
    diag(Sml) <- 1
    Sml <- diag(sds) %*% Sml %*% diag(sds)
    return(Sml)
  }
  if (!is.null(covMat)){
    # covariance known up to a constant
    c <- sum(diag(solve(covMat) %*% (t(Ys) %*% Ys)/nrow(Ys))) / ncol(Ys)
    Sml <- c * covMat
    return(Sml)
  }
  if (!is.null(corMat)){
    # correlation matrix known, variances unknown
    V <- eigen(corMat)
    D <- abs(V$values)
    V <- V$vectors
    Vinner <- eigen(diag(1/sqrt(D)) %*% t(V) %*% ((t(Ys) %*% Ys)/nrow(Ys))
                    %*% V %*% diag(1/sqrt(D)))
    Vinner$values <- sqrt(abs(Re(Vinner$values)))
    sds <- Re(diag(V %*% diag(sqrt(D)) %*% Vinner$vectors
                   %*% diag(Vinner$values) %*% t(Vinner$vectors)
                   %*% diag(sqrt(D)) %*% t(V)))
    Sml <- diag(sds) %*% corMat %*% diag(sds)
    return(Sml)
  }
  if (!is.null(covMat) & !is.null(corMat) & corType=="none" & varType=="none"){
    return(covML(Ys))
  }
}









#' Evaluate numerical properties square matrix
#' 
#' Function that evaluates various numerical properties of a square input
#' matrix. The intended use is to evaluate the various numerical properties of
#' what is assumed to be a covariance matrix. Another use is to evaluate the
#' various numerical properties of a (regularized) precision matrix.
#' 
#' The function evaluates various numerical properties of a covariance or
#' precision input matrix. The function assesses if the input matrix is
#' symmetric, if all its eigenvalues are real, if all its eigenvalues are
#' strictly positive, and if it is a diagonally dominant matrix. In addition,
#' the function calculates the trace, the determinant, and the spectral
#' condition number of the input matrix. See, e.g., Harville (1997) for more
#' details on the mentioned (numerical) matrix properties.
#' 
#' @param S Covariance or (regularized) precision \code{matrix}.
#' @param verbose A \code{logical} indicating if output should be printed on
#' screen.
#' @return \item{symm}{A \code{logical} indicating if the matrix is symmetric.}
#' \item{realEigen}{A \code{logical} indicating if the eigenvalues are real.}
#' \item{posEigen}{A \code{logical} indicating if the eigenvalues are strictly
#' positive.} \item{dd}{A \code{logical} indicating if the matrix is diagonally
#' dominant.} \item{trace}{A \code{numerical} giving the value of the trace.}
#' \item{det}{A \code{numerical} giving the value of the determinant.}
#' \item{condNumber}{A \code{numerical} giving the value of the spectral
#' condition number.}
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}}
#' @references Harville, D.A.(1997). Matrix algebra from a statistician's
#' perspective. New York: Springer-Verlag.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Evaluate numerical properties covariance matrix
#' ## Obtain, e.g., value trace
#' Seval <- evaluateS(Cx); Seval
#' Seval$trace
#' 
#' ## Evaluate numerical properties precision matrix after regularization
#' P <- ridgeP(Cx, lambda = 10, type = 'Alt')
#' Peval <- evaluateS(P); Peval
#' 
#' @export evaluateS
evaluateS <- function(S, verbose = TRUE){
  ##############################################################################
  # - Function evualuating various properties of an input matrix
  # - Intended use is to evaluate the various properties of what is assumed to
  #   be a covariance matrix
  # - Another use is to evaluate the various properties of a (regularized)
  #   precision matrix
  # - S       > sample covariance/correlation matrix or (regularized) precision
  #             matrix
  # - verbose > logical indicating if output should be printed on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(S)){
    stop("S should be a matrix")
  }
  else if (nrow(S) != ncol(S)){
    stop("S should be a square matrix")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    Sproperties <- list()

    # Is S symmetric?
    Sproperties$symm <- isSymmetric(S)

    # Are eigenvalues S real and positive?
    evs                   <- eigen(S)$values
    Sproperties$realEigen <- all(Im(evs) == 0)
    Sproperties$posEigen  <- all(evs > 0)

    # Is S diagonally dominant?
    Sproperties$diagDom <- all(abs(cov2cor(S)[upper.tri(S)]) < 1)

    # Trace and determinant S
    Sproperties$trace <- sum(diag(S))
    Sproperties$det   <- det(S)

    # Spectral condition number S
    Sproperties$condNumber <- abs(max(evs) / min(evs))

    if (verbose){
      cat("Properties of input matrix:\n")
      cat("----------------------------------------\n")
      cat("       symmetric : ", Sproperties$symm, "\n", sep="")
      cat("eigenvalues real : ", Sproperties$realEigen, "\n", sep="")
      cat(" eigenvalues > 0 : ", Sproperties$posEigen, "\n", sep="")
      cat("  diag. dominant : ", Sproperties$diagDom, "\n\n", sep="")
      cat("           trace : ", round(Sproperties$trace, 5), "\n", sep="")
      cat("     determinant : ", round(Sproperties$det, 5), "\n", sep="")
      cat(" l2 cond. number : ", round(Sproperties$condNumber, 5), "\n", sep="")
      cat("----------------------------------------\n")
    }

    # Return
    return(Sproperties)
  }
}









#' Compute partial correlation matrix or standardized precision matrix
#' 
#' Function computing the partial correlation matrix or standardized precision
#' matrix from an input precision matrix.
#' 
#' The function assumes that the input \code{matrix} is a precision matrix. If
#' \code{pc = FALSE} the standardized precision matrix, rather than the partial
#' correlation matrix, is given as the output value. The standardized precision
#' matrix is equal to the partial correlation matrix up to the sign of
#' off-diagonal entries.
#' 
#' @param P (Possibly regularized) precision \code{matrix}.
#' @param pc A \code{logical} indicating if the partial correlation matrix
#' should be computed.
#' @return A partial correlation \code{matrix} or a standardized precision
#' \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain regularized precision matrix
#' P <- ridgeP(Cx, lambda = 10, type = "Alt")
#' 
#' ## Obtain partial correlation matrix
#' pcor(P)
#' 
#' @export pcor
pcor <- function(P, pc = TRUE){
  ##############################################################################
  # - Function computing partial correlation/standardized precision matrix from
  #   a precision matrix
  # - P  > precision matrix (possibly regularized inverse of covariance or
  #        correlation matrix)
  # - pc > logical indicating if the partial correlation matrix should be
  #        computed
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(P)){
    stop("P should be a matrix")
  }
  else if (!isSymmetric(P)){
    stop("P should be a symmetric matrix")
  }
  else {
    # Compute partial correlation matrix
    if (pc){
      P       <- -P
      diag(P) <- -diag(P)
      Pcor    <- cov2cor(P)
      return(Pcor)
    }

    # Compute standardized precision matrix
    else {
      SP <- cov2cor(P)
      return(SP)
    }
  }
}









#' Generate a (data-driven) default target for usage in ridge-type shrinkage
#' estimation
#' 
#' Function that generates a (data-driven) default target for usage in (type I)
#' ridge shrinkage estimation of the precision matrix (see
#' \code{\link{ridgeP}}). The target that is generated is to be understood in
#' precision terms. Most options for target generation result in a target that
#' implies a situation of rotation equivariant estimation (see
#' \code{\link{ridgeP}}).
#' 
#' The function can generate the following default target matrices: \itemize{
#' \item \code{DAIE}: Diagonal matrix with average of inverse nonzero
#' eigenvalues of S as entries; \item \code{DIAES}: Diagonal matrix with
#' inverse of average of eigenvalues of S as entries; \item \code{DUPV}:
#' Diagonal matrix with unit partial variance as entries (identity matrix);
#' \item \code{DAPV}: Diagonal matrix with average of inverse variances of
#' \code{S} as entries; \item \code{DCPV}: Diagonal matrix with constant
#' partial variance as entries. Allows one to use other constant than DAIE,
#' DIAES, DUPV, DAPV, and in a sense Null; \item \code{DEPV}: Diagonal matrix
#' with the inverse variances of \code{S} as entries; \item \code{Null}: Null
#' matrix. } The targets \code{DUPV}, \code{DCPV}, and \code{Null} are not
#' data-driven in the sense that the input matrix \code{S} only provides
#' information on the size of the desired target. The targets \code{DAIE},
#' \code{DIAES}, \code{DAPV}, and \code{DEPV} are data-driven in the sense that
#' the input matrix \code{S} provides the information for the diagonal entries.
#' The argument \code{fraction} is only used when \code{type = "DAIE"}. The
#' argument \code{const} is only used when \code{type = "DCPV"}. All types
#' except \code{DEPV} and \code{Null} lead to rotation equivariant alternative
#' and archetypal Type I ridge estimators. The target \code{Null} also leads to
#' a rotation equivariant alternative Type II ridge estimator (see
#' \code{\link{ridgeP}}). Note that the \code{DIAES}, \code{DAPV}, and
#' \code{DEPV} targets amount to the identity matrix when the sample covariance
#' matrix \code{S} is standardized to be the correlation matrix. The same goes,
#' naturally, for the \code{DCPV} target when \code{const} is specified to be
#' 1.
#' 
#' @param S Sample covariance \code{matrix}.
#' @param type A \code{character} determining the type of default target. Must
#' be one of: "DAIE", "DIAES", "DUPV", "DAPV", "DCPV", "DEPV", "Null".
#' @param fraction A \code{numeric} indicating the fraction of the largest
#' eigenvalue below which an eigenvalue is considered zero.
#' @param const A \code{numeric} constant representing the partial variance.
#' @return Function returns a target \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}}
#' @references van Wieringen, W.N. & Peeters, C.F.W. (2016).  Ridge Estimation
#' of Inverse Covariance Matrices from High-Dimensional Data, Computational
#' Statistics & Data Analysis, vol. 103: 284-303.  Also available as
#' arXiv:1403.0904v3 [stat.ME].
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain default diagonal target matrix
#' default.target(Cx)
#' 
#' @export default.target
default.target <- function(S, type = "DAIE", fraction = 1e-04, const){
  ##############################################################################
  # - Function that generates a (data-driven) default target for usage in
  #   ridge-type shrinkage estimation
  # - The target that is generated is to be understood in precision terms
  # - See function 'ridgeS'
  # S        > sample covariance/correlation matrix
  # type     > character determining the type of default target;
  #            default = "DAIE" (see notes below)
  # fraction > fraction of largest eigenvalue below which an eigenvalue is
  #            considered zero. Only when type = "DAIE".
  # const    > numeric constant that represents the partial variance.
  #            Only when type = "DCPV"
  #
  # Notes:
  # - The van Wieringen-Peeters type I estimator and the archetypal I estimator
  #   utilize a p.d. target
  # - DAIE: diagonal average inverse eigenvalue
  #   Diagonal matrix with average of inverse nonzero eigenvalues of S as
  #   entries
  # - DIAES: diagonal inverse average eigenvalue S
  #   Diagonal matrix with inverse of average of eigenvalues of S as entries
  # - DUPV: diagonal unit partial variance
  #   Diagonal matrix with unit partial variance as entries (identity matrix)
  # - DAPV: diagonal average partial variance
  #   Diagonal matrix with average of inverse variances of S as entries
  # - DCPV: diagonal constant partial variance
  #   Diagonal matrix with constant partial variance as entries. Allows one to
  #   use other constant than [DAIE, DUPV, DAPV, and in a sense Null]
  # - DEPV: diagonal empirical partial variance
  #   Diagonal matrix with the inverse of variances of S as entries
  # - Null: Null matrix
  #   Matrix with only zero entries
  # - All but DEPV and Null lead to rotation equivariant alternative and
  #   archetype I ridge estimators
  # - Null also leads to a rotation equivariant alternative Type II estimator
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (class(type) != "character"){
    stop("Input (type) is of wrong class")
  }
  else if (!(type %in% c("DAIE", "DIAES", "DUPV", "DAPV",
                         "DCPV", "DEPV", "Null"))){
    stop("Input (type) should be one of {'DAIE', 'DIAES', 'DUPV', 'DAPV', ",
         "'DCPV', 'DEPV', 'Null'}")
  }
  else {
    # Compute and return a default target matrix
    # Diagonal matrix with average of inverse nonzero eigenvalues of S as
    # entries
    if (type == "DAIE"){
      if (class(fraction) != "numeric"){
        stop("Input (fraction) is of wrong class")
      } else if (length(fraction) != 1){
        stop("Length input (fraction) must be one")
      } else if (fraction < 0 | fraction > 1){
        stop("Input (fraction) is expected to be in the interval [0,1]")
      } else {
        Eigs   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        const  <- mean(1/(Eigs[Eigs >= Eigs[1]*fraction]))
        target <- const * diag(ncol(S))
      }
    }

    # Diagonal matrix with inverse of average of eigenvalues of S as entries
    if (type == "DIAES"){
      Eigs   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      const  <- 1/mean(Eigs)
      target <- const * diag(ncol(S))
    }

    # Diagonal matrix with unit partial variance as entries
    if (type == "DUPV"){
      target <- diag(ncol(S))
    }

    # Diagonal matrix with average empirical partial variances as entries
    if (type == "DAPV"){
      apv    <- mean(1/diag(S))
      target <- apv * diag(ncol(S))
    }

    # Diagonal matrix with constant partial variance as entries
    if (type == "DCPV"){
      if (class(const) != "numeric"){
        stop("Input (const) is of wrong class")
      } else if (length(const) != 1){
        stop("Length input (const) must be one")
      } else if (const <= 0 | const > .Machine$double.xmax){
        stop("Input (const) is expected to be in the interval (0, Inf)")
      } else {
        target <- const * diag(ncol(S))
      }
    }

    # Diagonal matrix with empirical partial variances as entries
    if (type == "DEPV"){
      target <- diag(1/diag(S))
    }

    # Null matrix
    if (type == "Null"){
      target <- matrix(0, ncol(S), nrow(S))
    }

    # Return
    colnames(target) = rownames(target) <- rownames(S)
    return(target)
  }
}




##------------------------------------------------------------------------------
##
## Function for Ridge Estimation of the Precision matrix
##
##------------------------------------------------------------------------------







#' Ridge estimation for high-dimensional precision matrices
#' 
#' Function that calculates various Ridge estimators for high-dimensional
#' precision matrices.
#' 
#' The function can calculate various ridge estimators for high-dimensional
#' precision matrices. Current (well-known) ridge estimators can be roughly
#' divided in two archetypes. The first archetypal form employs a convex
#' combination of \eqn{\mathbf{S}} and a positive definite (p.d.) target matrix
#' \eqn{\mathbf{T}}:
#' \eqn{\hat{\mathbf{\Omega}}^{\mathrm{I}}(\lambda_{\mathrm{I}}) =
#' [(1-\lambda_{\mathrm{I}}) \mathbf{S} + \lambda_{\mathrm{I}}
#' \mathbf{T}]^{-1}}, with \eqn{\lambda_{\mathrm{I}} \in (0,1]}. A common
#' target choice is for \eqn{\mathbf{T}} to be diagonal with
#' \eqn{(\mathbf{T})_{jj} = (\mathbf{S})_{jj}} for \eqn{j=1, \ldots, p}. The
#' second archetypal form can be given as
#' \eqn{\hat{\mathbf{\Omega}}^{\mathrm{II}}(\lambda_{\mathrm{II}}) =
#' (\mathbf{S} + \lambda_{\mathrm{II}} \mathbf{I}_{p})^{-1}} with
#' \eqn{\lambda_{\mathrm{II}} \in (0, \infty)}. Viewed from a penalized
#' estimation perspective, the two archetypes utilize penalties that do not
#' coincide with the matrix-analogue of the common ridge penalty. van Wieringen
#' and Peeters (2015) derive analytic expressions for alternative Type I and
#' Type II ridge precision estimators based on a proper L2-penalty. Their
#' alternative Type I estimator (target shrinkage) takes the form
#' \deqn{\hat{\mathbf{\Omega}}^{\mathrm{I}a}(\lambda_{a}) =
#' \left\{\left[\lambda_{a}\mathbf{I}_{p} + \frac{1}{4}(\mathbf{S} -
#' \lambda_{a}\mathbf{T})^{2}\right]^{1/2} + \frac{1}{2}(\mathbf{S} -
#' \lambda_{a}\mathbf{T})\right\}^{-1},} while their alternative Type II
#' estimator can be given as a special case of the former:
#' \deqn{\hat{\mathbf{\Omega}}^{\mathrm{II}a}(\lambda_{a}) =
#' \left\{\left[\lambda_{a}\mathbf{I}_{p} +
#' \frac{1}{4}\mathbf{S}^{2}\right]^{1/2} +
#' \frac{1}{2}\mathbf{S}\right\}^{-1}.} These alternative estimators were shown
#' to be superior to the archetypes in terms of risk under various loss
#' functions (van Wieringen and Peeters, 2015).
#' 
#' The \code{lambda} parameter in \code{ridgeP} generically indicates the
#' penalty parameter. It must be chosen in accordance with the type of ridge
#' estimator employed. The domains for the penalty parameter in the archetypal
#' estimators are given above. The domain for \code{lambda} in the alternative
#' estimators is \eqn{(0, \infty)}. The \code{type} parameter specifies the
#' type of ridge estimator. Specifying \code{type = "ArchI"} leads to usage of
#' the archetypal I estimator while specifying \code{type = "ArchII"} leads to
#' usage of the archetypal II estimator. In the latter situation the argument
#' \code{target} remains unused. Specifying \code{type = "Alt"} enables usage
#' of the alternative ridge estimators: when \code{type = "Alt"} and the
#' \code{target} matrix is p.d. one obtains the alternative Type I estimator;
#' when \code{type = "Alt"} and the \code{target} matrix is specified to be the
#' null-matrix one obtains the alternative Type II estimator.
#' 
#' The Type I estimators thus employ target shrinkage. The default target for
#' both the archetype and alternative is \code{default.target(S)}. When
#' \code{target} is not the null-matrix it is expected to be p.d. for the
#' alternative type I estimator. The target is always expected to be p.d. in
#' case of the archetypal I estimator. The archetypal Type I ridge estimator is
#' rotation equivariant when the target is of the form \eqn{\mu\mathbf{I}_{p}}
#' with \eqn{\mu \in (0,\infty)}. The archetypal Type II estimator is rotation
#' equivariant by definition. When the target is of the form
#' \eqn{\varphi\mathbf{I}_{p}} with \eqn{\varphi \in [0,\infty)}, then the
#' alternative ridge estimator is rotation equivariant. Its analytic
#' computation is then particularly speedy as the (relatively) expensive matrix
#' square root can then be circumvented.
#' 
#' @param S Sample covariance \code{matrix}.
#' @param lambda A \code{numeric} representing the value of the penalty
#' parameter.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @return Function returns a regularized precision \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Anders E. Bilgrau
#' @seealso \code{\link{default.target}}
#' @references van Wieringen, W.N. & Peeters, C.F.W. (2016). Ridge Estimation
#' of Inverse Covariance Matrices from High-Dimensional Data, Computational
#' Statistics & Data Analysis, vol. 103: 284-303. Also available as
#' arXiv:1403.0904v3 [stat.ME].
#' 
#' van Wieringen, W.N. & Peeters, C.F.W. (2015). Application of a New Ridge
#' Estimator of the Inverse Covariance Matrix to the Reconstruction of
#' Gene-Gene Interaction Networks. In: di Serio, C., Lio, P., Nonis, A., and
#' Tagliaferri, R. (Eds.) `Computational Intelligence Methods for
#' Bioinformatics and Biostatistics'. Lecture Notes in Computer Science, vol.
#' 8623. Springer, pp. 170-179.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain regularized precision matrix
#' ridgeP(Cx, lambda = 10, type = "Alt")
#' 
#' @export ridgeP
ridgeP <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Function that calculates ridge estimators of a precision matrix
  # - Version that uses the RcppArmadillo implementation
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a precision
  #             matrix
  # - ArchI   > Archetypal I ridge estimator of a precision matrix
  # - ArchII  > Archetypal II ridge estimator of a precision matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II estimator
  # - When target is not the null-matrix it is expected to be p.d. for the
  #   vWP type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  if (!isSymmetric(S)) {
    stop("S should be a symmetric matrix")
  }
  else if (lambda <= 0) {
    stop("lambda should be positive")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else{
    # Calculate Ridge estimator
    # Alternative estimator
    if (type == "Alt"){
      if (!isSymmetric(target)) {
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]) {
        stop("S and target should be of the same dimension")
      } else {
        P_Alt <- .armaRidgeP(S, target, lambda)
      }
      dimnames(P_Alt) <- dimnames(S)
      return(symm(P_Alt))
    }

    # Archetypal I
    if (type == "ArchI"){
      if (lambda > 1){
        stop("lambda should be in (0,1] for this type of Ridge estimator")
      } else if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (any(eigen(target, symmetric = TRUE,
                           only.values = TRUE)$values <= 0)){
        stop("Target should always be p.d. for this type of ridge estimator")
      } else {
        P_ArchI <- solve((1-lambda) * S + lambda * solve(target))
        return(symm(P_ArchI))
      }
    }

    # Archetypal II
    if (type == "ArchII"){
      P_ArchII <- solve(S + lambda * diag(nrow(S)))
      return(symm(P_ArchII))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Penalty Parameter selection
##
##------------------------------------------------------------------------------







#' Select optimal penalty parameter by approximate leave-one-out
#' cross-validation
#' 
#' Function that selects the optimal penalty parameter for the
#' \code{\link{ridgeP}} call by usage of approximate leave-one-out
#' cross-validation. Its output includes (a.o.) the precision matrix under the
#' optimal value of the penalty parameter.
#' 
#' The function calculates an approximate leave-one-out cross-validated
#' (aLOOCV) negative log-likelihood score (using a regularized ridge estimator
#' for the precision matrix) for each value of the penalty parameter contained
#' in the search grid. The utilized aLOOCV score was proposed by Lian (2011)
#' and Vujacic et al. (2014). The aLOOCV negative log-likeliho od score is
#' computationally more efficient than its non-approximate counterpart (see
#' \code{\link{optPenalty.LOOCV}}). For details on the aLOOCV negative
#' log-likelihood score see Lian (2011) and Vujacic et al (2014). For scalar
#' matrix targets (see \code{\link{default.target}}) the complete solution path
#' of the alternative Type I and II ridge estimators (see \code{\link{ridgeP}})
#' depends on only 1 eigendecomposition and 1 matrix inversion, making the
#' determination of the optimal penalty value particularly efficient (see van
#' Wieringen and Peeters, 2015).
#' 
#' The value of the penalty parameter that achieves the lowest aLOOCV negative
#' log-likelihood score is deemed optimal. The penalty parameter must be
#' positive such that \code{lambdaMin} must be a positive scalar. The maximum
#' allowable value of \code{lambdaMax} depends on the type of ridge estimator
#' employed. For details on the type of ridge estimator one may use (one of:
#' "Alt", "ArchI", "ArchII") see \code{\link{ridgeP}}. The ouput consists of an
#' object of class list (see below). When \code{output = "light"} (default)
#' only the \code{optLambda} and \code{optPrec} elements of the list are given.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param step An \code{integer} determining the number of steps in moving
#' through the grid [\code{lambdaMin}, \code{lambdaMax}].
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param cor A \code{logical} indicating if the evaluation of the approximate
#' LOOCV score should be performed on the correlation scale.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param output A \code{character} indicating if the output is either heavy or
#' light. Must be one of: "all", "light".
#' @param graph A \code{logical} indicating if the grid search for the optimal
#' penalty parameter should be visualized.
#' @param verbose A \code{logical} indicating if information on progress should
#' be printed on screen.
#' @return An object of class list: \item{optLambda}{A \code{numeric} giving
#' the optimal value of the penalty parameter.} \item{optPrec}{A \code{matrix}
#' representing the precision matrix of the chosen type (see
#' \code{\link{ridgeS}}) under the optimal value of the penalty parameter.}
#' \item{lambdas}{A \code{numeric} vector representing all values of the
#' penalty parameter for which approximate cross-validation was performed; Only
#' given when \code{output = "all"}.} \item{aLOOCVs}{A \code{numeric} vector
#' representing the approximate cross-validated negative log-likelihoods for
#' each value of the penalty parameter given in \code{lambdas}; Only given when
#' \code{output = "all"}.}
#' @note When \code{cor = TRUE} correlation matrices are used in the
#' computation of the approximate (cross-validated) negative log-likelihood
#' score, i.e., the sample covariance matrix is a matrix on the correlation
#' scale. When performing evaluation on the correlation scale the data are
#' assumed to be standardized. If \code{cor = TRUE} and one wishes to used the
#' default target specification one may consider using \code{target =
#' default.target(covML(Y, cor = TRUE))}. This gives a default target under the
#' assumption of standardized data.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.LOOCV}},
#' \code{\link{optPenalty.LOOCVauto}}, \cr \code{\link{default.target}},
#' \code{\link{covML}}
#' @references Lian, H. (2011). Shrinkage tuning parameter selection in
#' precision matrices estimation. Journal of Statistical Planning and
#' Inference, 141: 2839-2848.
#' 
#' van Wieringen, W.N. & Peeters, C.F.W. (2016). Ridge Estimation of Inverse
#' Covariance Matrices from High-Dimensional Data, Computational Statistics &
#' Data Analysis, vol. 103: 284-303. Also available as arXiv:1403.0904v3
#' [stat.ME].
#' 
#' Vujacic, I., Abbruzzo, A., and Wit, E.C. (2014). A computationally fast
#' alternative to cross-validation in penalized Gaussian graphical models.
#' arXiv: 1309.6216v2 [stat.ME].
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT  <- optPenalty.aLOOCV(X, lambdaMin = .001, lambdaMax = 30, step = 400); OPT
#' OPT$optLambda	# Optimal penalty
#' OPT$optPrec	  # Regularized precision under optimal penalty
#' 
#' ## Another example with standardized data
#' X <- scale(X, center = TRUE, scale = TRUE)
#' OPT  <- optPenalty.aLOOCV(X, lambdaMin = .001, lambdaMax = 30,
#'                           step = 400, cor = TRUE,
#'                           target = default.target(covML(X, cor = TRUE))); OPT
#' OPT$optLambda	# Optimal penalty
#' OPT$optPrec	  # Regularized precision under optimal penalty
#' 
#' @export optPenalty.aLOOCV
optPenalty.aLOOCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt",
                              cor = FALSE, target = default.target(covML(Y)),
                              output = "light", graph = TRUE, verbose = TRUE) {
  ##############################################################################
  # - Function that selects the optimal penalty parameter by approximate
  #   leave-one-out cross-validation
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax   > maximum value penalty parameter (dependent on 'type')
  # - step        > determines the coarseness in searching the grid
  #                 [lambdaMin, lambdaMax]
  # - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - cor         > logical indicating if evaluation of the LOOCV score should be
  #                 performed on the correlation matrix
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - output      > must be one of {"all", "light"}, default = "light"
  # - graph       > Optional argument for visualization optimal penalty
  #                 selection, default = TRUE
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")
  # require("sfsmisc")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (class(cor) != "logical"){
    stop("Input (cor) is of wrong class")
  }
  else if (!(output %in% c("all", "light"))){
    stop("Input (output) should be one of {'all', 'light'}")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else {
    # Set preliminaries
    S       <- covML(Y, cor = cor)
    n       <- nrow(Y)
    lambdas <- lseq(lambdaMin, lambdaMax, length = step)
    aLOOCVs <- numeric()

    # Calculate approximate LOOCV scores
    if (verbose){cat("Calculating approximate LOOCV scores...", "\n")}
    if (type == "Alt" & all(target == 0)){
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Covariance matrix based on data input (Y) and target should be ",
             "of the same dimension")
      } else {
        Spectral <- eigen(S, symmetric = TRUE)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k])
          P         <- Spectral$vectors %*% diag(1/Eigshrink) %*% t(Spectral$vectors)
          nLL       <- .5 * .LL(S, P)
          isum      <- numeric()

          for (i in 1:n){
            S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
            isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
          }

          aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
          if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
        }
      }
    } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0) &
               (length(unique(diag(target))) == 1)) {
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Covariance matrix based on data input (Y) and target should be ",
             "of the same dimension")
      } else if (any(diag(target) <= 0)){
        stop("Input (target) should be p.d.")
      } else {
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k], const = varPhi)
          P         <- Spectral$vectors %*% diag(1/Eigshrink) %*% t(Spectral$vectors)
          nLL       <- .5 * .LL(S, P)
          isum      <- numeric()

          for (i in 1:n){
            S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
            isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
          }

          aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
          if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
        }
      }
    } else {
      for (k in 1:length(lambdas)){
        P    <- ridgeP(S, lambdas[k], type = type, target = target)
        nLL  <- .5 * .LL(S, P)
        isum <- numeric()

        for (i in 1:n){
          S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
          isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
        }

        aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
        if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
      }
    }

    # Visualization
    optLambda <- min(lambdas[which(aLOOCVs == min(aLOOCVs))])
    if (graph){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
      plot(log(lambdas), type = "l", aLOOCVs, axes = FALSE,
           xlab = "ln(penalty value)",
           ylab = "Approximate LOOCV neg. log-likelihood", main = Main)
      axis(2, ylim = c(min(aLOOCVs),max(aLOOCVs)), col = "black", lwd = 1)
      axis(1, col = "black", lwd = 1)
      par(xpd = FALSE)
      abline(h = min(aLOOCVs), v = log(optLambda), col = "red")
      legend("topright",
             legend = c(paste("min. approx. LOOCV neg. LL: ",
                              round(min(aLOOCVs),3), sep = ""),
                        paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
    }

    # Return
    if (output == "all"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target),
                  lambdas = lambdas, aLOOCVs = aLOOCVs))
    }
    if (output == "light"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target)))
    }
  }
}









#' Select optimal penalty parameter by \eqn{K}-fold cross-validation
#' 
#' Function that selects the optimal penalty parameter for the
#' \code{\link{ridgeP}} call by usage of \eqn{K}-fold cross-validation. Its
#' output includes (a.o.) the precision matrix under the optimal value of the
#' penalty parameter.
#' 
#' The function calculates a cross-validated negative log-likelihood score
#' (using a regularized ridge estimator for the precision matrix) for each
#' value of the penalty parameter contained in the search grid by way of
#' \eqn{K}-fold cross-validation. The value of the penalty parameter that
#' achieves the lowest cross-validated negative log-likelihood score is deemed
#' optimal. The penalty parameter must be positive such that \code{lambdaMin}
#' must be a positive scalar. The maximum allowable value of \code{lambdaMax}
#' depends on the type of ridge estimator employed. For details on the type of
#' ridge estimator one may use (one of: "Alt", "ArchI", "ArchII") see
#' \code{\link{ridgeP}}. The ouput consists of an object of class list (see
#' below). When \code{output = "light"} (default) only the \code{optLambda} and
#' \code{optPrec} elements of the list are given.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param step An \code{integer} determining the number of steps in moving
#' through the grid [\code{lambdaMin}, \code{lambdaMax}].
#' @param fold A \code{numeric} or \code{integer} specifying the number of
#' folds to apply in the cross-validation.
#' @param cor A \code{logical} indicating if the evaluation of the LOOCV score
#' should be performed on the correlation scale.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param output A \code{character} indicating if the output is either heavy or
#' light. Must be one of: "all", "light".
#' @param graph A \code{logical} indicating if the grid search for the optimal
#' penalty parameter should be visualized.
#' @param verbose A \code{logical} indicating if information on progress should
#' be printed on screen.
#' @return An object of class list: \item{optLambda}{A \code{numeric} giving
#' the optimal value of the penalty parameter.} \item{optPrec}{A \code{matrix}
#' representing the precision matrix of the chosen type (see
#' \code{\link{ridgeP}}) under the optimal value of the penalty parameter.}
#' \item{lambdas}{A \code{numeric} vector representing all values of the
#' penalty parameter for which cross-validation was performed; Only given when
#' \code{output = "all"}.} \item{LLs}{A \code{numeric} vector representing the
#' mean of cross-validated negative log-likelihoods for each value of the
#' penalty parameter given in \code{lambdas}; Only given when \code{output =
#' "all"}.}
#' @note When \code{cor = TRUE} correlation matrices are used in the
#' computation of the (cross-validated) negative log-likelihood score, i.e.,
#' the \eqn{K}-fold sample covariance matrix is a matrix on the correlation
#' scale. When performing evaluation on the correlation scale the data are
#' assumed to be standardized. If \code{cor = TRUE} and one wishes to used the
#' default target specification one may consider using \code{target =
#' default.target(covML(Y, cor = TRUE))}. This gives a default target under the
#' assumption of standardized data.
#' 
#' Under the default setting of the fold-argument, \code{fold = nrow(Y)}, one
#' performes leave-one-out cross-validation.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.kCVauto}},
#' \code{\link{optPenalty.aLOOCV}}, \cr \code{\link{default.target}},
#' \code{\link{covML}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty using K = n
#' OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100); OPT
#' OPT$optLambda	# Optimal penalty
#' OPT$optPrec	  # Regularized precision under optimal penalty
#' 
#' ## Another example with standardized data
#' X <- scale(X, center = TRUE, scale = TRUE)
#' OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100, cor = TRUE,
#'                        target = default.target(covML(X, cor = TRUE))); OPT
#' OPT$optLambda	# Optimal penalty
#' OPT$optPrec	  # Regularized precision under optimal penalty
#' 
#' ## Another example using K = 5
#' OPT  <- optPenalty.kCV(X, lambdaMin = .5, lambdaMax = 30, step = 100, fold = 5); OPT
#' OPT$optLambda	# Optimal penalty
#' OPT$optPrec	  # Regularized precision under optimal penalty
#' 
#' @export optPenalty.kCV
optPenalty.kCV <- function(Y, lambdaMin, lambdaMax, step, fold = nrow(Y),
                           cor = FALSE, target = default.target(covML(Y)),
                           type = "Alt", output = "light", graph = TRUE,
                           verbose = TRUE) {
  ##############################################################################
  # - Function that selects the optimal penalty parameter by leave-one-out
  #   cross-validation
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax   > maximum value penalty parameter (dependent on 'type')
  # - step        > determines the coarseness in searching the grid
  #                 [lambdaMin, lambdaMax]
  # - fold        > cross-validation fold, default gives LOOCV
  # - cor         > logical indicating if evaluation of the LOOCV score should be
  #                 performed on the correlation matrix
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - output      > must be one of {"all", "light"}, default = "light"
  # - graph       > Optional argument for visualization optimal penalty
  #                 selection, default = TRUE
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("graphics")
  # require("sfsmisc")

  if (class(verbose) != "logical")
  { stop("Input (verbose) is of wrong class") }
  if (verbose){ cat("Perform input checks...", "\n") }
  if (!is.matrix(Y))
  { stop("Input (Y) should be a matrix") }
  if (class(lambdaMin) != "numeric")
  { stop("Input (lambdaMin) is of wrong class") }
  if (length(lambdaMin) != 1)
  { stop("Input (lambdaMin) must be a scalar") }
  if (lambdaMin <= 0)
  { stop("Input (lambdaMin) must be positive") }
  if (class(lambdaMax) != "numeric")
  { stop("Input (lambdaMax) is of wrong class") }
  if (length(lambdaMax) != 1)
  { stop("Input (lambdaMax) must be a scalar") }
  if (lambdaMax <= lambdaMin)
  { stop("Input (lambdaMax) must be larger than lambdaMin") }
  if (class(step) != "numeric")
  { stop("Input (step) is of wrong class") }
  if (!.is.int(step))
  { stop("Input (step) should be integer") }
  if (step <= 0)
  { stop("Input (step) should be a positive integer") }
  if (class(cor) != "logical")
  { stop("Input (cor) is of wrong class") }
  if (!(output %in% c("all", "light")))
  { stop("Input (output) should be one of {'all', 'light'}") }
  if (class(graph) != "logical")
  { stop("Input (graph) is of wrong class") }
  if (class(fold) != "numeric" & class(fold) != "integer")
  { stop("Input (fold) is of wrong class") }
  if ((fold <=  1) | (fold > nrow(Y)))
  { stop("Input (fold) out of range") }

  # make k-folds as list
  fold <- max(min(ceiling(fold), nrow(Y)), 2)
  fold <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
  shuffle <- sample(1:nrow(Y), nrow(Y))
  folds <- split(shuffle, as.factor(fold))

  # Set preliminaries
  LLs     <- numeric()
  lambdas <- lseq(lambdaMin, lambdaMax, length=step)

  # Calculate CV scores
  if (verbose) {
    cat("Calculating cross-validated negative log-likelihoods...\n")
  }
  for (k in 1:length(lambdas)){
    LLs <- c(LLs, .kcvl(lambdas[k], Y, cor, target, type, folds))
    if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
  }

  # Visualization
  optLambda <- min(lambdas[which(LLs == min(LLs))])
  if (graph){
    if (type == "Alt"){Main = "Alternative ridge estimator"}
    if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
    if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
    plot(log(lambdas), type = "l", LLs, axes = FALSE,
         xlab = "ln(penalty value)",
         ylab = "LOOCV neg. log-likelihood", main = Main)
    axis(2, ylim = c(min(LLs),max(LLs)), col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    par(xpd = FALSE)
    abline(h = min(LLs), v = log(optLambda), col = "red")
    legend("topright",
           legend = c(paste("min. LOOCV neg. LL: ", round(min(LLs),3),sep=""),
                      paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
  }

  # Return
  S <- covML(Y, cor = cor)
  if (output == "all"){
    return(list(optLambda = optLambda,
                optPrec = ridgeP(S, optLambda, type = type, target = target),
                lambdas = lambdas, LLs = LLs))
  }
  if (output == "light"){
    return(list(optLambda = optLambda,
                optPrec = ridgeP(S, optLambda, type = type, target = target)))
  }
}









#' Automatic search for optimal penalty parameter
#' 
#' Function that performs an 'automatic' search for the optimal penalty
#' parameter for the \code{\link{ridgeP}} call by employing Brent's method to
#' the calculation of a cross-validated negative log-likelihood score.
#' 
#' The function determines the optimal value of the penalty parameter by
#' application of the Brent algorithm (1971) to the \eqn{K}-fold
#' cross-validated negative log-likelihood score (using a regularized ridge
#' estimator for the precision matrix). The search for the optimal value is
#' automatic in the sense that in order to invoke the root-finding abilities of
#' the Brent method, only a minimum value and a maximum value for the penalty
#' parameter need to be specified as well as a starting penalty value. The
#' value at which the \eqn{K}-fold cross-validated negative log-likelihood
#' score is minimized is deemed optimal. The function employs the Brent
#' algorithm as implemented in the
#' \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.htmloptim}
#' function.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param lambdaInit A \code{numeric} giving the initial (starting) value for
#' the penalty parameter.
#' @param fold A \code{numeric} or \code{integer} specifying the number of
#' folds to apply in the cross-validation.
#' @param cor A \code{logical} indicating if the evaluation of the LOOCV score
#' should be performed on the correlation scale.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @return An object of class \code{list}: \item{optLambda}{A \code{numeric}
#' giving the optimal value for the penalty parameter.} \item{optPrec}{A
#' \code{matrix} representing the precision matrix of the chosen type (see
#' \code{\link{ridgeP}}) under the optimal value of the penalty parameter.}
#' @note When \code{cor = TRUE} correlation matrices are used in the
#' computation of the (cross-validated) negative log-likelihood score, i.e.,
#' the \eqn{K}-fold sample covariance matrix is a matrix on the correlation
#' scale. When performing evaluation on the correlation scale the data are
#' assumed to be standardized. If \code{cor = TRUE} and one wishes to used the
#' default target specification one may consider using \code{target =
#' default.target(covML(Y, cor = TRUE))}. This gives a default target under the
#' assumption of standardized data.
#' 
#' Under the default setting of the fold-argument, \code{fold = nrow(Y)}, one
#' performes leave-one-out cross-validation.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{GGMblockNullPenalty}}, \code{\link{GGMblockTest}},
#' \code{\link{ridgeP}}, \code{\link{optPenalty.aLOOCV}},
#' \code{\link{optPenalty.kCV}}, \cr \code{\link{default.target}},
#' \code{\link{covML}}
#' @references Brent, R.P. (1971). An Algorithm with Guaranteed Convergence for
#' Finding a Zero of a Function. Computer Journal 14: 422-425.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty using K = n
#' OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30); OPT
#' OPT$optLambda # Optimal penalty
#' OPT$optPrec   # Regularized precision under optimal penalty
#' 
#' ## Another example with standardized data
#' X <- scale(X, center = TRUE, scale = TRUE)
#' OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30, cor = TRUE,
#'                           target = default.target(covML(X, cor = TRUE))); OPT
#' OPT$optLambda # Optimal penalty
#' OPT$optPrec   # Regularized precision under optimal penalty
#' 
#' ## Another example using K = 5
#' OPT <- optPenalty.kCVauto(X, lambdaMin = .001, lambdaMax = 30, fold = 5); OPT
#' OPT$optLambda # Optimal penalty
#' OPT$optPrec   # Regularized precision under optimal penalty
#' 
#' @export optPenalty.kCVauto
optPenalty.kCVauto <- function(Y, lambdaMin, lambdaMax,
                               lambdaInit = (lambdaMin + lambdaMax)/2,
                               fold = nrow(Y), cor = FALSE,
                               target = default.target(covML(Y)),
                               type = "Alt") {
    ##############################################################################
    # - Function that determines the optimal value of the penalty parameter by
    #   application of the Brent algorithm to the (leave-one-out) cross-validated
    #   log-likelihood
    # - Y          > (raw) Data matrix, variables in columns
    # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
    # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
    # - lambdaInit > initial value for lambda for starting optimization
    # - fold       > cross-validation fold, default gives LOOCV
    # - cor        > logical indicating if evaluation of the LOOCV score should be
    #                performed on the correlation matrix
    # - target     > target (precision terms) for Type I estimators,
    #                default = default.target(covML(Y))
    # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
    ##############################################################################

    # Dependencies
    # require("base")
    # require("stats")

    # input checks
    if (!is.matrix(Y))
      { stop("Input (Y) should be a matrix") }
    if (class(lambdaMin) != "numeric")
      { stop("Input (lambdaMin) is of wrong class") }
    if (length(lambdaMin) != 1)
      { stop("Input (lambdaMin) must be a scalar") }
    if (lambdaMin <= 0)
      { stop("Input (lambdaMin) must be positive") }
    if (class(lambdaMax) != "numeric")
      { stop("Input (lambdaMax) is of wrong class") }
    if (length(lambdaMax) != 1)
      { stop("Input (lambdaMax) must be a scalar") }
    if (lambdaMax <= lambdaMin)
      { stop("Input (lambdaMax) must be larger than lambdaMin") }
    if (class(lambdaInit) != "numeric")
      { stop("Input (lambdaInit) is of wrong class") }
    if (length(lambdaInit) != 1)
      { stop("Input (lambdaInit) must be a scalar") }
    if (lambdaInit <= lambdaMin)
      { stop("Input (lambdaInit) must be larger than lambdaMin") }
    if (lambdaInit > lambdaMax)
      { stop("Input (lambdaInit) must be smaller than lambdaMax") }
    if (class(cor) != "logical")
      { stop("Input (cor) is of wrong class") }
    if (class(fold) != "numeric" & class(fold) != "integer")
      { stop("Input (fold) is of wrong class") }
    if ((fold <=  1) | (fold > nrow(Y)))
      { stop("Input (fold) out of range") }

    # make k-folds as list
    fold <- max(min(ceiling(fold), nrow(Y)), 2)
    fold <- rep(1:fold, ceiling(nrow(Y)/fold))[1:nrow(Y)]
    shuffle <- sample(1:nrow(Y), nrow(Y))
    folds <- split(shuffle, as.factor(fold))

    # determine optimal value of ridge penalty parameter
    optLambda <- optim(lambdaInit, .kcvl, method = "Brent", lower = lambdaMin,
                       upper = lambdaMax, Y = Y, cor = cor, target = target,
                       type = type, folds = folds)$par

    # Return
    return(list(optLambda = optLambda,
                optPrec = ridgeP(covML(Y, cor = cor), optLambda,
                                 type = type, target = target)))
}









#' Visualize the spectral condition number against the regularization parameter
#' 
#' Function that visualizes the spectral condition number of the regularized
#' precision matrix against the domain of the regularization parameter. The
#' function can be used to heuristically determine an acceptable (minimal)
#' value for the penalty parameter.
#' 
#' Under certain target choices the proposed ridge estimators (see
#' \code{\link{ridgeP}}) are rotation equivariant, i.e., the eigenvectors of
#' \eqn{\mathbf{S}} are left intact. Such rotation equivariant situations help
#' to understand the effect of the ridge penalty on the precision estimate: The
#' effect can be understood in terms of shrinkage of the eigenvalues of the
#' unpenalized precision estimate \eqn{\mathbf{S}^{-1}}. Maximum shrinkage
#' implies that all eigenvalues are forced to be equal (in the rotation
#' equivariant situation). The spectral condition number w.r.t. inversion
#' (ratio of maximum to minimum eigenvalue) of the regularized precision matrix
#' may function as a heuristic in determining the (minimal) value of the
#' penalty parameter. A matrix with a high condition number is near-singular
#' (the relative distance to the set of singular matrices equals the reciprocal
#' of the condition number; Demmel, 1987) and its inversion is numerically
#' unstable. Such a matrix is said to be ill-conditioned. Numerically,
#' ill-conditioning will mean that small changes in the penalty parameter lead
#' to dramatic changes in the condition number. From a numerical point of view
#' one can thus track the domain of the penalty parameter for which the
#' regularized precision matrix is ill-conditioned. When plotting the condition
#' number against the (domain of the) penalty parameter, there is a point of
#' relative stabilization (when working in the \eqn{p > n} situation) that can
#' be characterized by a leveling-off of the acceleration along the curve when
#' plotting the condition number against the (chosen) domain of the penalty
#' parameter. This suggest the following fast heuristic for determining the
#' (minimal) value of the penalty parameter: The value of the penalty parameter
#' for which the spectral condition number starts to stabilize may be termed an
#' acceptable (minimal) value.
#' 
#' The function outputs a graph of the (spectral) matrix condition number over
#' the domain [\code{lambdaMin}, \code{lambdaMax}]. When \code{norm = "2"} the
#' spectral condition number is calculated. It is determined by exact
#' calculation using the spectral decomposition. For most purposes this exact
#' calculation is fast enough, especially when considering rotation equivariant
#' situations (see \code{\link{ridgeP}}). For such situations the amenities for
#' fast eigenvalue calculation as provided by
#' \href{https://CRAN.R-project.org/package=RSpectraRSpectra} are used
#' internally. When exact computation of the spectral condition number is
#' deemed too costly one may approximate the computationally friendly
#' L1-condition number. This approximation is accessed through the
#' \link[base:kappa]{rcond} function (Anderson et al. 1999).
#' 
#' When \code{Iaids = TRUE} the basic condition number plot is amended/enhanced
#' with two additional plots (over the same domain of the penalty parameter as
#' the basic plot): The approximate loss in digits of accuracy (for the
#' operation of inversion) and an approximation to the second-order derivative
#' of the curvature in the basic plot. These interpretational aids can enhance
#' interpretation of the basic condition number plot and may support choosing a
#' value for the penalty parameter (see Peeters, van de Wiel, & van Wieringen,
#' 2016). When \code{vertical = TRUE} a vertical line is added at the constant
#' \code{value}. This option can be used to assess if the optimal penalty
#' obtained by, e.g., the routines \code{\link{optPenalty.LOOCV}} or
#' \code{\link{optPenalty.aLOOCV}}, has led to a precision estimate that is
#' well-conditioned.
#' 
#' @param S Sample covariance \code{matrix}.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param step An \code{integer} determining the number of steps in moving
#' through the grid [\code{lambdaMin}, \code{lambdaMax}].
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param norm A \code{character} indicating the norm under which the condition
#' number is to be calculated/estimated. Must be one of: "1", "2".
#' @param Iaids A \code{logical} indicating if the basic condition number plot
#' should be amended with interpretational aids.
#' @param vertical A \code{logical} indicating if output graph should come with
#' a vertical line at a pre-specified value for the penalty parameter.
#' @param value A \code{numeric} indicating a pre-specified value for the
#' penalty parameter.
#' @param main A \code{character} with which to specify the main title of the
#' output graph.
#' @param nOutput A \code{logical} indicating if numeric output should be
#' returned.
#' @param verbose A \code{logical} indicating if information on progress should
#' be printed on screen.
#' @param suppressChecks A \code{logical} indicating if the input checks should
#' be suppressed.
#' @return The function returns a graph. If \code{nOutput = TRUE} the function
#' also returns an object of class \code{list}: \item{lambdas}{A \code{numeric}
#' vector representing all values of the penalty parameter for which the
#' condition number was calculated. The values of the penalty parameter are
#' log-equidistant.} \item{conditionNumbers}{A \code{numeric} vector containing
#' the condition number for each value of the penalty parameter given in
#' \code{lambdas}.}
#' @note The condition number of a (regularized) covariance matrix is
#' equivalent to the condition number of its corresponding inverse, the
#' (regularized) precision matrix. Please note that the \code{target} argument
#' (for Type I ridge estimators) is assumed to be specified in precision terms.
#' In case one is interested in the condition number of a Type I estimator
#' under a covariance target, say \eqn{\mathbf{\Gamma}}, then just specify
#' \code{target = solve}(\eqn{\mathbf{\Gamma}}).
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}},
#' \code{\link{optPenalty.LOOCV}}, \code{\link{optPenalty.aLOOCV}},
#' \code{\link{default.target}}
#' @references Anderson, E, Bai, Z., ..., Sorenson, D. (1999). LAPACK Users'
#' Guide (3rd ed.). Philadelphia, PA, USA: Society for Industrial and Applied
#' Mathematics.
#' 
#' Demmel, J.W. (1987). On condition numbers and the distance to the nearest
#' ill-posed problem. Numerische Mathematik, 51: 251--289.
#' 
#' Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2020). The
#' spectral condition number plot for regularization parameter evaluation.
#' Computational Statistics, 35: 629--646.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Assess spectral condition number across grid of penalty parameter
#' CNplot(Cx, lambdaMin = .0001, lambdaMax = 50, step = 1000)
#' 
#' ## Include interpretational aids
#' CNplot(Cx, lambdaMin = .0001, lambdaMax = 50, step = 1000, Iaids = TRUE)
#' 
#' @export CNplot
CNplot <- function(S, lambdaMin, lambdaMax, step, type = "Alt",
                   target = default.target(S, type = "DUPV"), norm = "2",
                   Iaids = FALSE, vertical = FALSE, value = 1e-100,
                   main = "", nOutput = FALSE, verbose = TRUE,
                   suppressChecks = FALSE){
  #############################################################################
  # - Function that visualizes the spectral condition number against the
  #   regularization parameter
  # - Can be used to heuristically determine the (minimal) value of the
  #   penalty parameter
  # - The ridge estimators operate by shrinking the eigenvalues
  # - This is especially the case when targets are used that lead to
  #   rotation equivariant estimators
  # - Maximum shrinkage (under rotation equivariance) implies that all
  #   eigenvalues will be equal
  # - Ratio of maximum and minimum eigenvalue of P can then function
  #   as a heuristic
  # - It's point of stabilization can give an acceptable value for the penalty
  # - The ratio boils down to the (spectral) condition number of a matrix
  # - S         > sample covariance/correlation matrix
  # - lambdaMin > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax > maximum value penalty parameter (dependent on 'type')
  # - step      > determines the coarseness in searching the grid
  #               [lambdaMin, lambdaMax]. The steps on the grid are equidistant
  #               on the log scale
  # - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target    > target (precision terms) for Type I estimators,
  #               default = default.target(S)
  # - norm      > indicates the norm under which the condition number is to be
  #               estimated. Default is the L2-norm. The L1-norm can be (cheaply)
  #               approximated
  # - Iaids     > logical indicating if interpretational aids should also be
  #               visualized. The aids are the approximate loss in digits of
  #               accuracy and an approximation of the acceleration along the
  #               curve. Default = FALSE
  # - vertical  > optional argument for visualization vertical line in graph
  #               output, default = FALSE. Can be used to indicate the value of,
  #               e.g., the optimal penalty as indicated by some routine. Can
  #               be used to assess if this optimal penalty will lead to a
  #               well-conditioned estimate
  # - value     > indicates constant on which to base vertical line when
  #               vertical = TRUE. Default is a very small value
  # - main      > logical indicating if plot should contain type of estimator
  #               as main title
  # - nOutput   > logical indicating if numeric output should be given
  #               (lambdas and condition numbers)
  # - verbose   > logical indicating if process information should be printed
  #               on-screen
  # - suppressChecks > logical indicating if input checks should be
  #               suppressed. When the variable-dimension is very large, the
  #               full set of input checks (supplied as a courtesy to the user)
  #               can take up considerable time. Suppressing these checks
  #               will speed up the procedure (when the user knows what he or
  #               she is doing).
  #############################################################################

  # Dependencies
  # require("base")
  # require("graphics")
  # require("Hmisc")
  # require("sfsmisc")
  # require("RSpectra")

  if (!suppressChecks){
    if (class(verbose) != "logical"){
      stop("Input (verbose) is of wrong class")
    }
    if (verbose){
      cat("Perform input checks...", "\n")
    }
    if (!is.matrix(S)){
      stop("S should be a matrix")
    }
    else if (class(lambdaMin) != "numeric"){
      stop("Input (lambdaMin) is of wrong class")
    }
    else if (length(lambdaMin) != 1){
      stop("lambdaMin must be a scalar")
    }
    else if (lambdaMin <= 0){
      stop("lambdaMin must be positive")
    }
    else if (class(lambdaMax) != "numeric"){
      stop("Input (lambdaMax) is of wrong class")
    }
    else if (length(lambdaMax) != 1){
      stop("lambdaMax must be a scalar")
    }
    else if (lambdaMax <= lambdaMin){
      stop("lambdaMax must be larger than lambdaMin")
    }
    else if (class(step) != "numeric"){
      stop("Input (step) is of wrong class")
    }
    else if (!.is.int(step)){
      stop("step should be integer")
    }
    else if (step <= 0){
      stop("step should be a positive integer")
    }
    else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
      stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
    }
    else if (dim(target)[1] != dim(S)[1]){
      stop("S and target should be of the same dimension")
    }
    else if (type == "ArchI" & lambdaMax > 1){
      stop("lambda should be in (0,1] for this type of Ridge estimator")
    }
    else if (!(norm %in% c("2", "1"))){
      stop("norm should be one of {'2', '1'}")
    }
    else if (class(Iaids) != "logical"){
      stop("Input (Iaids) is of wrong class")
    }
    else if (class(vertical) != "logical"){
      stop("Input (vertical) is of wrong class")
    }
    else if (vertical & (class(value) != "numeric")){
      stop("Input (value) is of wrong class")
    }
    else if (vertical & (any(value <= 0))){
      stop("Input (value) must be strictly positive")
    }
    else if (vertical & (length(value) != 1)){
      stop("Input (value) must be a scalar")
    }
    else if (class(main) != "character"){
      stop("Input (main) is of wrong class")
    }
    else if (class(nOutput) != "logical"){
      stop("Input (nOutput) is of wrong class")
    }
  }

  # Set preliminaries
  lambdas <- lseq(lambdaMin, lambdaMax, length = step)
  condNR  <- numeric()

  if (norm == "2"){
    # Calculate spectral condition number ridge estimate on lambda grid
    if (verbose){cat("Calculating spectral condition numbers...", "\n")}
    if (type == "Alt" & all(target == 0)){
      Spectral <- eigs_sym(S, 2, which = "BE",
                           opts = list(retvec = FALSE, maxitr = 1000000))$values
      for (k in 1:length(lambdas)){
        Eigshrink <- .armaEigShrink(Spectral, lambdas[k])
        condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
      }
    } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0) &
               (length(unique(diag(target))) == 1)){
      varPhi   <- unique(diag(target))
      Spectral <- eigs_sym(S, 2, which = "BE",
                           opts = list(retvec = FALSE, maxitr = 1000000))$values
      for (k in 1:length(lambdas)){
        Eigshrink <- .armaEigShrink(Spectral, lambdas[k], cons = varPhi)
        condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
      }
    } else {
      if (type == "Alt"){
        for (k in 1:length(lambdas)){
          Eigs      <- .armaEigShrinkAnyTarget(S, target = target, lambdas[k])
          condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
        }
      } else if (type == "ArchI" & all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
        varPhi   <- unique(diag(target))
        Spectral <- eigs_sym(S, 2, which = "BE",
                             opts = list(retvec = FALSE, maxitr = 1000000))$values
        for (k in 1:length(lambdas)){
          Eigshrink <- .armaEigShrinkArchI(Spectral, lambdas[k], cons = varPhi)
          condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
        }
      } else {
        if (type == "ArchI"){
          for (k in 1:length(lambdas)){
            P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
            Eigs      <- eigs_sym(P, 2, which = "BE",
                                  opts = list(retvec = FALSE, maxitr = 1000000))$values
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
        if (type == "ArchII"){
          Spectral <- eigs_sym(S, 2, which = "BE",
                               opts = list(retvec = FALSE, maxitr = 1000000))$values
          for (k in 1:length(lambdas)){
            Eigs      <- Spectral + lambdas[k]
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
      }
    }
  }

  if (norm == "1"){
    # Calculate approximation to condition number under 1-norm
    if (verbose){cat("Approximating condition number under 1-norm...", "\n")}
    if (type == "Alt"){
      for (k in 1:length(lambdas)){
        P         <- .armaRidgeP(S, target = target, lambdas[k])
        condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
      }
    }
    if (type != "Alt"){
      for (k in 1:length(lambdas)){
        P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
        condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
      }
    }
  }

  if (Iaids) {
    # Make calculations for interpretational aids
    if (verbose){cat("Calculating interpretational aids...", "\n")}
    dLoss   <- floor(log10(condNR))
    logLamb <- log(lambdas)
    delta   <- logLamb[2] - logLamb[1]
    which   <- c(1, length(condNR))
    Core    <- condNR[-which]
    up      <- condNR[-c(which[2]-1, which[2])]
    down    <- condNR[-c(which[1], which[1]+1)]
    cdapp   <- (down - (2 * Core) + up)/(delta^2)
  }

  # Visualization
  if (verbose){cat("Plotting...", "\n")}
  if (norm == "2"){Ylab = "spectral condition number"}
  if (norm == "1"){Ylab = "condition number under 1-norm"}
  if (Iaids){par(mfrow=c(1,3))}
  plot(log(lambdas), type = "l", condNR, axes = FALSE, col = "blue4",
       xlab = "ln(penalty value)", ylab = Ylab, main = main)
  axis(2, ylim = c(0,max(condNR)), col = "black", lwd = 1)
  axis(1, col = "black", lwd = 1)
  minor.tick(nx = 10, ny = 0, tick.ratio = .4)
  if (vertical){abline(v = log(value), col = "red")}
  par(xpd = FALSE)
  if (Iaids){
    plot(log(lambdas), dLoss, axes = FALSE, type = "l",
         col = "green3", xlab = "ln(penalty value)",
         ylab = "approximate loss in digits of accuracy")
    axis(2, ylim = c(0,max(dLoss)), col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    minor.tick(nx = 10, ny = 0, tick.ratio = .4)
    if (vertical){abline(v = log(value), col = "red")}
    xlimits <- range(log(lambdas))
    plot(log(lambdas[-which]), cdapp, axes = FALSE, type = "l",
         col = "orange", xlim = xlimits, xlab = "ln(penalty value)",
         ylab = "approximation of acceleration")
    axis(2, ylim = c(0,max(cdapp)), col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    minor.tick(nx = 10, ny = 0, tick.ratio = .4)
    if (vertical){abline(v = log(value), col = "red")}
    par(mfrow=c(1,1))
  }

  # Possible output
  if (nOutput){
    return(list(lambdas = lambdas, conditionNumbers = condNR))
  }
}




##------------------------------------------------------------------------------
##
## Functions for Block Independence Testing and Mutual Information
##
##------------------------------------------------------------------------------

if(getRversion() >= "2.15.1") utils::globalVariables("rags2ridges")







#' Generate the distribution of the penalty parameter under the null hypothesis
#' of block-independence
#' 
#' Function that serves as a precursor function to the block-independence test
#' (see \code{\link{GGMblockTest}}). It generates an empirical distribution of
#' the penalty parameter under the null hypothesis of block independence (in
#' the regularized precision matrix).
#' 
#' This function can be viewed as a precursor to the function for the
#' block-independence test (see \code{\link{GGMblockTest}}). The mentioned test
#' evaluates the null hypothesis of block-independence against the alternative
#' of block-dependence (presence of non-zero elements in the off-diagonal
#' block) in the precision matrix using high-dimensional data. To accommodate
#' the high-dimensionality the parameters of interest are estimated in a
#' penalized manner (ridge-type penalization, see \code{\link{ridgeP}}).
#' Penalization involves a degree of freedom (the penalty parameter) which
#' needs to be fixed before testing. This function then generates an empirical
#' distribution of this penalty parameter. Hereto the samples are permutated
#' within block. The resulting permuted data sets represent the null
#' hypothesis. To avoid the dependence on a single permutation, many permuted
#' data sets are generated. For each permutation the optimal penalty parameter
#' is determined by means of cross-validation (see
#' \code{\link{optPenalty.LOOCVauto}}). The resulting optimal penalty
#' parameters are returned. An estimate of the location (such as the median) is
#' recommended for use in the block-independence test.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param id A \code{numeric} vector acting as an indicator variable for two
#' blocks of the precision matrix. The blocks should be coded as \code{0} and
#' \code{1}.
#' @param nPerm A \code{numeric} or \code{integer} determining the number of
#' permutations.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param lambdaInit A \code{numeric} giving the initial value for the penalty
#' parameter for starting optimization.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param ncpus A \code{numeric} or \code{integer} indicating the desired
#' number of cpus to be used.
#' @return A \code{numeric} vector, representing the distribution of the (LOOCV
#' optimal) penalty parameter under the null hypothesis of block-independence.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.LOOCVauto}},
#' \code{\link{default.target}}, \code{\link{GGMblockTest}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 15
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:15] = letters[1:15]
#' id <- c(rep(0, 10), rep(1, 5))
#' 
#' ## Generate null distribution of the penalty parameter
#' lambda0dist <- GGMblockNullPenalty(X, id, 5, 0.001, 10)
#' 
#' ## Location of null distribution
#' lambdaNull <- median(lambda0dist)
#' 
#' @export GGMblockNullPenalty
GGMblockNullPenalty <- function(Y, id, nPerm = 25, lambdaMin, lambdaMax,
                                lambdaInit = (lambdaMin+lambdaMax)/2,
                                target = default.target(covML(Y)),
                                type = "Alt", ncpus = 1){
  ##############################################################################
  # - Function generating the distribution of the penalty parameter
  # - It does so under the null hypothesis of block independence
  # - The optimal value of the penalty parameter is determined under multiple
  #   permutations
  # - Optimal penalty determined using the 'optPenalty.LOOCVauto' function
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - nPerm      > desired number of permutations
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - ncpus      > desired number of cpus
  #
  # Notes:
  # - Dependency on 'snowfall' when ncpus > 1
  ##############################################################################

  # Dependencies
  # require("base")
  # require("snowfall")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input (Y) should not contain missings")
  }
  else if (class(id) != "numeric" & class(id) != "integer"){
    stop("Input (id) is of wrong class")
  }
  else if (!(all(unique(id) %in% c(0, 1)))){
    stop("Input (id) has unlawful entries")
  }
  else if (length(id) != ncol(Y)){
    stop("Column dimension input (Y) does not match with length input (id)")
  }
  else if (class(nPerm) != "numeric" & class(nPerm) != "integer"){
    stop("Input (nPerm) is of wrong class")
  }
  else if (!.is.int(nPerm)){
    stop("Input (nPerm) is expected to be a (numeric) integer")
  }
  else if (nPerm <= 0){
    stop("Input (nPerm) must be strictly positive")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be strictly positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than input (lambdaMin)")
  }
  else if (class(lambdaInit) != "numeric"){
    stop("Input (lambdaInit) is of wrong class")
  }
  else if (length(lambdaInit) != 1){
    stop("Input (lambdaInit) must be a scalar")
  }
  else if (lambdaInit <= lambdaMin){
    stop("Input (lambdaInit) must be larger than input (lambdaMin)")
  }
  else if (lambdaMax <= lambdaInit){
    stop("Input (lambdaInit) must be smaller than input (lambdaMax)")
  }
  else if (class(ncpus) != "numeric" & class(ncpus) != "integer"){
    stop("Input (ncpus) is of wrong class")
  }
  else if (!.is.int(ncpus)){
    stop("Input (ncpus) is expected to be a (numeric) integer")
  }
  else if (ncpus <= 0){
    stop("Input (ncpus) must be strictly positive")
  }
  else {
    # Determine null distribution
    if (ncpus == 1){
      lambdaNull <- sapply(1:nPerm, .lambdaNullDist, Y = Y, id = id,
                           lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                           lambdaInit = lambdaInit, target = target,
                           type = type)
    }
    if (ncpus > 1){
      sfInit(parallel = TRUE, cpus = ncpus)
      sfLibrary(rags2ridges, verbose = FALSE)
      lambdaNull <- sfSapply(1:nPerm, .lambdaNullDist, Y = Y, id = id,
                             lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                             lambdaInit = lambdaInit, target = target,
                             type = type)
      sfStop()
    }

    # Return
    return(lambdaNull)
  }
}









#' Test for block-indepedence
#' 
#' Function performing a test that evaluates the null hypothesis of
#' block-independence against the alternative of block-dependence (presence of
#' non-zero elements in the off-diagonal block) in the precision matrix using
#' high-dimensional data. The mentioned test is a permutation-based test (see
#' details).
#' 
#' The function performs a permutation test for the null hypothesis of
#' block-independence against the alternative of block-dependence (presence of
#' non-zero elements in the off-diagonal block) in the precision matrix using
#' high-dimensional data. In the low-dimensional setting the common test
#' statistic under multivariate normality (cf. Anderson, 2003) is: \deqn{ \log(
#' \| \hat{\mathbf{\Sigma}}_a \| ) + \log( \| \hat{\mathbf{\Sigma}}_b \| ) -
#' \log( \| \hat{\mathbf{\Sigma}} \| ), } where the
#' \eqn{\hat{\mathbf{\Sigma}}_a}, \eqn{\hat{\mathbf{\Sigma}}_b},
#' \eqn{\hat{\mathbf{\Sigma}}} are the estimates of the covariance matrix in
#' the sub- and whole group(s), respectively.
#' 
#' To accommodate the high-dimensionality the parameters of interest are
#' estimated in a penalized manner (ridge-type penalization, see
#' \code{\link{ridgeP}}). Penalization involves a degree of freedom (the
#' penalty parameter: \code{lambda}) which needs to be fixed before testing. To
#' decide on the penalty parameter for testing we refer to the
#' \code{\link{GGMblockNullPenalty}} function. With an informed choice on the
#' penalty parameter at hand, the null hypothesis is evaluated by permutation.
#' Hereto the samples are permutated within block. The resulting permuted data
#' set represents the null hypothesis. Many permuted data sets are generated.
#' For each permutation the test statistic is calculated. The observed test
#' statistic is compared to the null distribution from the permutations.
#' 
#' The function implements an efficient permutation resampling algorithm (see
#' van Wieringen et al., 2008, for details.): If the probability of a p-value
#' being below \code{lowCiThres} is smaller than 0.001 (read: the test is
#' unlikely to become significant), the permutation analysis is terminated and
#' a p-value of unity (1) is reported.
#' 
#' When \code{verbose = TRUE} also graphical output is generated: A histogram
#' of the null-distribution. Note that, when \code{ncpus} is larger than 1,
#' functionalities from
#' \href{https://cran.r-project.org/package=snowfallsnowfall} are imported.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param id A \code{numeric} vector acting as an indicator variable for two
#' blocks of the precision matrix. The blocks should be coded as \code{0} and
#' \code{1}.
#' @param nPerm A \code{numeric} or \code{integer} determining the number of
#' permutations.
#' @param lambda A \code{numeric} representing the penalty parameter employed
#' in the permutation test.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param lowCiThres A \code{numeric} taking a value between 0 and 1.
#' Determines speed of efficient p-value calculation.
#' @param ncpus A \code{numeric} or \code{integer} indicating the desired
#' number of cpus to be used.
#' @param verbose A \code{logical} indicating if information on progress and
#' output should be printed on screen.
#' @return An object of class list: \item{statistic}{A \code{numeric}
#' representing the observed test statistic (i.e., likelihood ratio).}
#' \item{pvalue}{A \code{numeric} giving the p-value for the block-independence
#' test.} \item{nulldist}{A \code{numeric} vector representing the permutation
#' null distribution for the test statistic.} \item{nperm}{A \code{numeric}
#' indicating the number of permutations used for p-value calculation.}
#' \item{remark}{A \code{"character"} that states whether the permutation
#' algorithm was terminated prematurely or not.}
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.LOOCVauto}},
#' \code{\link{default.target}}, \code{\link{GGMblockNullPenalty}}
#' @references Anderson, T.W. (2003). An Introduction to Multivariate
#' Statistical Analysis, 3rd Edition. John Wiley.
#' 
#' van Wieringen, W.N., van de Wiel, M.A., and van der Vaart, A.W. (2008). A
#' Test for Partial Differential Expression. Journal of the American
#' Statistical Association 103: 1039-1049.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 15
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:15] = letters[1:15]
#' id <- c(rep(0, 10), rep(1, 5))
#' 
#' ## Generate null distribution of the penalty parameter
#' lambda0dist <- GGMblockNullPenalty(X, id, 5, 0.001, 10)
#' 
#' ## Location of null distribution
#' lambdaNull <- median(lambda0dist)
#' 
#' ## Perform test
#' testRes <- GGMblockTest(X, id, nPerm = 100, lambdaNull)
#' 
#' @export GGMblockTest
GGMblockTest <- function (Y, id, nPerm = 1000, lambda,
                          target = default.target(covML(Y)), type = "Alt",
                          lowCiThres = 0.1, ncpus = 1, verbose = TRUE) {
  ##############################################################################
  # - Function performing a permutation test for block structure in the
  #   precision matrix
  # - The setting is a high-dimensional one
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - nPerm      > desired number of permutations, default = 1000
  # - lambda     > the penalty parameter employed in the permutation test
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - lowCiThres > A value between 0 and 1, determining speed of efficient
  #                p-value calculation
  # - ncpus      > desired number of cpus
  # - verbose    > logical indicating if progress/output should be printed on
  #                screen
  #
  # Notes:
  # - Dependency on 'snowfall' when ncpus > 1
  # - When verbose = TRUE, also graphical output is given: a histogram of the
  #   null-distribution
  # - The value for 'lambda' ideally stems from 'GGMblockNullPenalty'
  # - If the probability of a p-value being below 'lowCiThres' is smaller than
  #   0.001 (meaning: the test is unlikely to become significant), the
  #   permutation analysis is terminated and a p-value of unity (1) is reported
  ##############################################################################

  # Dependencies
  # require("base")
  # require("snowfall")
  # require("graphics")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input (Y) should not contain missings")
  }
  else if (class(id) != "numeric" & class(id) != "integer"){
    stop("Input (id) is of wrong class")
  }
  else if (!(all(unique(id) %in% c(0, 1)))){
    stop("Input (id) has unlawful entries")
  }
  else if (length(id) != ncol(Y)){
    stop("Column dimension input (Y) does not match with length input (id)")
  }
  else if (class(nPerm) != "numeric" & class(nPerm) != "integer"){
    stop("Input (nPerm) is of wrong class")
  }
  else if (!.is.int(nPerm)){
    stop("Input (nPerm) is expected to be a (numeric) integer")
  }
  else if (nPerm <= 0){
    stop("Input (nPerm) must be strictly positive")
  }
  else if (class(lambda) != "numeric"){
    stop("Input (lambda) is of wrong class")
  }
  else if (length(lambda) != 1){
    stop("Input (lambda) must be a scalar")
  }
  else if (lambda <= 0){
    stop("Input (lambda) must be strictly positive")
  }
  else if (class(lowCiThres) != "numeric"){
    stop("Input (lowCiThres) is of wrong class")
  }
  else if (length(lowCiThres) != 1){
    stop("Input (lowCiThres) must be a scalar")
  }
  else if (lowCiThres <= 0 | lowCiThres >= 1){
    stop("Input (lowCiThres) must be in the interval (0,1)")
  }
  else if (class(ncpus) != "numeric" & class(ncpus) != "integer"){
    stop("Input (ncpus) is of wrong class")
  }
  else if (!.is.int(ncpus)){
    stop("Input (ncpus) is expected to be a (numeric) integer")
  }
  else if (ncpus <= 0){
    stop("Input (ncpus) must be strictly positive")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # Observed test statistics
    S     <- solve(ridgeP(covML(Y), lambda = lambda,
                          target = target, type = type))
    llObs <- log(det(S[id == 0, id == 0, drop = FALSE])) +
      log(det(S[id == 1, id == 1, drop = FALSE])) - log(det(S))

    # Define steps at which the possibility of significance should be evaluated
    steps <- sort(unique(c(0, 25, 50, 100, 150, 200,
                           seq(from = 250, to = 2750, by = 250),
                           seq(from = 3000, to = 10000, by = 500),
                           seq(from = 11000, to = 50000, by = 1000), nPerm)))
    steps <- steps[steps <= nPerm]

    # Generate null distribution
    nullDist <- numeric()
    if (ceiling(ncpus) > 1){
      sfInit(parallel = TRUE, cpus = ncpus)
      sfLibrary(rags2ridges, verbose = FALSE)
    }

    for (j in 1:(length(steps) - 1)){
      if (verbose){
        cat(paste(steps[j], " of ", steps[length(steps)],
                  " permutations done, and counting...\n", sep = ""))
      }
      if (ncpus == 1){
        nullDistPart <- sapply(c((steps[j] + 1):steps[j + 1]), .blockTestStat,
                               Y = Y, id = id, lambda = lambda, target = target,
                               type = type)
      }
      if (ncpus > 1){
        nullDistPart <- sfSapply(c((steps[j] + 1):steps[j + 1]), .blockTestStat,
                                 Y = Y, id = id, lambda = lambda,
                                 target = target, type = type)
      }
      nullDist <- c(nullDist, nullDistPart); rm(nullDistPart); gc()
      pVal     <- sum(nullDist >= as.numeric(llObs))/steps[j + 1]
      pBound   <- pVal - sqrt(pVal * (1 - pVal)/steps[j + 1]) * 3.09
      significanceUnlikely <- (pBound > lowCiThres)
      if (significanceUnlikely){pVal <- 1; break}
      if (verbose){cat(paste(steps[j + 1], "of", steps[length(steps)],
                             " permutations done", sep = " "), "\n")}
    }

    if (ncpus > 1){sfStop()}

    # Generating on screen (graphical) output
    if (verbose){
      # Visual summary of test results
      xlims     <- c(min(c(llObs, nullDist)), max(c(llObs, nullDist)))
      histFreqs <- hist(nullDist, n = sqrt(sum(nPerm))+1, plot = FALSE)$counts
      hist(nullDist, xlim = xlims, n = sqrt(sum(nPerm))+1, col = "blue",
           border = "lightblue", xlab = "null statistics",
           ylab = "frequency", main = "Histogram of null distribution")
      lines(rep(llObs, max(histFreqs)), 0.9 * (c(1:max(histFreqs))-1),
            col = "red", lwd = 2)
      text(quantile(c(nullDist, llObs), probs = 0.05), 0.95 * max(histFreqs),
           paste("p-value:", pVal))
      text(llObs, 0.95 * max(histFreqs), "test stat.")

      # Summary of test results
      remark <-
        ifelse(significanceUnlikely,
               "resampling terminated prematurely due to unlikely significance",
               "none")
      cat("\n")
      cat("Likelihood ratio test for block independence\n")
      cat("----------------------------------------\n")
      cat("-> number of permutations : ", nPerm, "\n", sep="")
      cat("-> test statistic         : ", round(llObs, digits = 3), "\n",sep="")
      cat("-> p-value                : ", round(pVal, digits = 3), "\n", sep="")
      cat("-> remark                 : ", remark, "\n", sep="")
      cat("----------------------------------------\n")
      cat("\n")
    }

    # Return
    return(list(statistic = llObs, pvalue = pVal, nulldist = nullDist,
                nperm = nPerm, remark = remark))
  }
}









#' Mutual information between two sets of variates within a multivariate normal
#' distribution
#' 
#' Function computing the mutual information between two exhaustive and
#' mutually exclusive splits of a set of multivariate normal random variables.
#' 
#' 
#' @param S A positive-definite covariance \code{matrix}.
#' @param split1 A \code{numeric}, indicating the variates (by column number)
#' forming the first split. The second split is automatically formed from its
#' complement.
#' @return A \code{numeric}, the mutual information between the variates
#' forming \code{split1} and those forming its complement.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}}.
#' @references Cover, T.M., Thomas, J.A. (2012), Elements of information
#' theory.
#' @examples
#' 
#' # create a covariance matrix
#' Sigma <- covML(matrix(rnorm(100), ncol=5))
#' 
#' # impulse response analysis
#' GGMmutualInfo(Sigma, c(1,2))
#' 
#' @export GGMmutualInfo
GGMmutualInfo <- function(S, split1){
  ##############################################################################
  # - Function that calculates the mutual information between two exhaustive and
  #   mutually exclusive splits of normal p-variate random variable.
  # - S      > Sample Covariance matrix.
  # - split1 > A numeric indicating the variates (by column number) forming the
  #            first split. The second split is automatically formed from its
  #            complement.
  #
  # NOTES:
  # - No dependencies at current
  ##############################################################################

  if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (length(split1) < 1 & length(split1) > nrow(S)-1){
    stop("Input (split1) is of wrong length.")
  }
  else {
    # mutual information
    MI <- log(det(S[-split1,-split1])) -
          log(det(S[-split1,-split1] - S[-split1,split1,drop=FALSE] %*%
                  solve(S[split1,split1,drop=FALSE])
                  %*% S[split1,-split1,drop=FALSE]))
    return(MI)
  }
}




##------------------------------------------------------------------------------
##
## Test for Vanishing Partial Correlations
##
##------------------------------------------------------------------------------







#' Determine the support of a partial correlation/precision matrix
#' 
#' Function that determines the support of a partial correlation/precision
#' matrix by thresholding and sparsifies it accordingly.
#' 
#' The function transforms the possibly regularized input precision matrix to a
#' partial correlation matrix. Subsequently, the support of this partial
#' correlation matrix is determined. Support determination is performed either
#' by simple thresholding on the absolute values of matrix entries
#' (\code{threshold = "absValue"}) or by usage of local FDR (\code{threshold =
#' "localFDR"}). A third option is to retain a prespecified number of matrix
#' entries based on absolute values. For example, one could wish to retain
#' those entries representing the ten strongest absolute partial correlations
#' (\code{threshold = "top"}). As a variation on this theme, a fourth option
#' (\code{threshold = "connected"}) retains the top edges such that the
#' resulting graph is connected (this may result in dense graphs in practice).
#' The argument \code{absValueCut} is only used when \code{threshold =
#' "absValue"}. The argument \code{top} is only used when \code{threshold =
#' "top"}. The argument \code{FDRcut} is only used when \code{threshold =
#' "localFDR"}.
#' 
#' The function is to some extent a wrapper around certain
#' \href{https://cran.r-project.org/package=fdrtoolfdrtool} functions when
#' \code{threshold = "localFDR"}. In that case a mixture model is fitted to the
#' nonredundant partial correlations by
#' \href{https://cran.r-project.org/package=fdrtoolfdrtool}. The decision to
#' retain elements is then based on the argument \code{FDRcut}. Elements with a
#' posterior probability \eqn{\geq} FDRcut (equalling 1 - local FDR) are
#' retained. See Schaefer and Strimmer (2005) for further details on usage of
#' local FDR in graphical modeling.
#' 
#' @param P (Possibly regularized) precision \code{matrix}.
#' @param threshold A \code{character} signifying type of sparsification by
#' thresholding. Must be one of: "absValue", "connected", "localFDR", "top".
#' @param absValueCut A \code{numeric} giving the cut-off for partial
#' correlation element selection based on absolute value thresholding.
#' @param FDRcut A \code{numeric} giving the cut-off for partial correlation
#' element selection based on local false discovery rate (FDR) thresholding.
#' @param top A \code{numeric} specifying the exact number of partial
#' correlation elements to retain based on absolute value.
#' @param output A \code{character} specifying the type of output required.
#' Must be one of: "heavy", "light".
#' @param verbose A \code{logical} indicating if intermediate output should be
#' printed on screen.
#' @return If the input \code{P} is a standardized precision (or partial
#' correlation) matrix the function returns a sparsified precision (or partial
#' correlation) \code{matrix} whenever \code{output = "heavy"}. If the input
#' \code{P} is an unstandardized precision matrix the function returns an
#' object of class \code{list} whenever \code{output = "heavy"}:
#' \item{sparseParCor}{A \code{matrix} representing the sparsified partial
#' correlation matrix.} \item{sparsePrecision}{A \code{matrix} representing the
#' sparsified precision matrix.}
#' 
#' When \code{output = "light"}, only the (matrix) positions of the zero and
#' non-zero elements are returned in an object of class \code{list}:
#' \item{zeros}{A \code{matrix} representing the row and column positions of
#' zero entries.} \item{nonzeros}{A \code{matrix} representing the row and
#' column positions of non-zero entries.}
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.aLOOCV}},
#' \code{\link{optPenalty.LOOCV}}
#' @references Schaefer, J., and Strimmer, K. (2005). A shrinkage approach to
#' large-scale covariance estimation and implications for functional genomics.
#' Statistical Applications in Genetics and Molecular Biology, 4:32.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)
#' 
#' ## Determine support regularized (standardized) precision under optimal penalty
#' sparsify(OPT$optPrec, threshold = "localFDR")
#' 
#' @export sparsify
sparsify <- function(P, threshold = c("absValue", "connected", "localFDR", "top"),
                     absValueCut = .25, FDRcut = .9, top = 10,
                     output = "heavy", verbose = TRUE){
  ##############################################################################
  # - Function that sparsifies/determines support of a partial correlation
  #   matrix
  # - Support can be determined by absolute value thresholding or by local FDR
  #   thresholding
  # - Local FDR operates on the nonredundant non-diagonal elements of a partial
  #   correlation matrix
  # - One can also choose to threshold based on the top X of absolute partial
  #   correlations
  # - One can also choose to threshold based on the minimum absolute partial
  #   correlation for which the resulting graph is connected
  # - Function is to some extent a wrapper around certain 'fdrtool' functions
  # - P           > (possibly shrunken) precision matrix
  # - threshold   > signifies type of thresholding
  # - absValueCut > cut-off for partial correlation elements selection based on
  #                 absolute value thresholding.
  #                 Only when threshold = 'absValue'. Default = .25
  # - FDRcut      > cut-off for partial correlation element selection based on
  #                 local FDR thresholding
  #                 Only when threshold = 'localFDR'. Default = .9
  # - top         > partial correlation element selection based on retainment
  #                 'top' number of absolute partial correlation.
  #                 Only when threshold = 'top'. Default = 10
  # - output      > must be one of {"heavy", "light"}, default = "heavy"
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen. Only when threshold = 'localFDR'. Default = TRUE.
  #
  # NOTES:
  # - Input (P) may be the partial correlation matrix or the standardized
  #   precision matrix. These are identical up to the signs of off-diagonal
  #   elements. Either can be used as it has no effect on the thresholding
  #   operator and the ensuing sparsified result.
  # - Input (P) may also be the unstandardized precision matrix. The function
  #   converts it to the partial correlation matrix
  # - The function evaluates if the input (P) is a partial
  #   correlation/standardized precision matrix or an unstandardized precision
  #   matrix. If the input amounts to the latter both the sparsified partial
  #   correlation matrix and the corresponding sparsified precision matrix are
  #   given as output (when output = "heavy"). Otherwise, the ouput consists of
  #   the sparsified partial correlation/standardized precision matrix.
  # - When output = "light", only the (matrix) positions of the zero and
  #   non-zero elements are returned.
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("fdrtool")
  # require("igraph")

  if (!is.matrix(P)){
    stop("Input (P) should be a matrix")
  }
  else if (!isSymmetric(P)){
    stop("Input (P) should be a symmetric matrix")
  }
  else if (!evaluateS(P, verbose = FALSE)$posEigen){
    stop("Input (P) is expected to be positive definite")
  }
  else if (missing(threshold)){
    stop("Need to specify type of sparsification ('absValue' or 'localFDR' ",
         "or 'connected' or 'top')")
  }
  else if (!(threshold %in% c("absValue", "connected", "localFDR", "top"))){
    stop("Input (threshold) should be one of
         {'absValue', 'connected', 'localFDR', 'top'}")
  }
  else if (!(output %in% c("light", "heavy"))){
    stop("Input (output) should be one of {'light', 'heavy'}")
  }
  else {
    # Obtain partial correlation matrix
    if (all(length(unique(diag(P))) == 1 & unique(diag(P)) == 1)){
      stan = TRUE
      PC  <- P
    } else {
      stan = FALSE
      PC  <- symm(pcor(P))
    }

    # Number of nonredundant elements
    NR <- (ncol(P)*(ncol(P) - 1))/2

    # Obtain sparsified matrix
    if (threshold == "top"){
      if (class(top) != "numeric"){
        stop("Input (top) is of wrong class")
      } else if (length(top) != 1){
        stop("Input (top) must be a scalar")
      } else if (!.is.int(top)){
        stop("Input (top) should be a numeric integer")
      } else if (top <= 0){
        stop("Input (top) must be strictly positive")
      } else if (top >= NR){
        stop("Input (top) must be smaller than the number of nonredundant ",
             "off-diagonal elements of the input matrix P")
      } else {
        absValueCut <- sort(abs(PC[upper.tri(PC)]),
                            decreasing = TRUE)[ceiling(top)]
        threshold   <- "absValue"
      }
    }

    if (threshold == "connected"){
      sumPC <- summary(abs(PC[upper.tri(PC)]))
      maxPC <- as.numeric(sumPC[6]); minPC <- as.numeric(sumPC[1])
      for (j in 1:100){
        absValueCut <- (maxPC + minPC)/2
        PC0 <- PC
        PC0[!(abs(PC0) >= absValueCut)] <- 0
        if (igraph::is.connected(graph.adjacency(adjacentMat(PC0), "undirected"))){
          minPC <- absValueCut
        } else {
          maxPC <- absValueCut
        }
        if (abs(absValueCut - (maxPC + minPC)/2) < 10^(-10)){
          absValueCut <- minPC; break
        }
      }
      threshold   <- "absValue"
    }

    if (threshold == "absValue"){
      if (class(absValueCut) != "numeric"){
        stop("Input (absValueCut) is of wrong class")
      } else if (length(absValueCut) != 1){
        stop("Input (absValueCut) must be a scalar")
      } else if (absValueCut <= 0 | absValueCut >= 1){
        stop("Input (absValueCut) must be in the interval (0,1)")
      } else {
        PC0 <- PC
        PC0[!(abs(PC0) >= absValueCut)] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    if (threshold == "localFDR"){
      if (class(FDRcut) != "numeric"){
        stop("Input (FDRcut) is of wrong class")
      } else if (length(FDRcut) != 1){
        stop("Input (FDRcut) must be a scalar")
      } else if (FDRcut <= 0 | FDRcut >= 1){
        stop("Input (FDRcut) must be in the interval (0,1)")
      } else if (class(verbose) != "logical"){
        stop("Input (verbose) is of wrong class")
      } else {
        lFDRs <- 1 - fdrtool(PC[upper.tri(PC)], "correlation",
                             plot = verbose, verbose = verbose)$lfdr
        PC0   <- diag(nrow(PC))
        PC0[lower.tri(PC0)] <- 1
        zeros <- which(PC0 == 0, arr.ind = TRUE)
        zeros <- zeros[which(lFDRs <= FDRcut),]
        PC0   <- PC
        PC0[zeros] <- 0
        PC0[cbind(zeros[,2], zeros[,1])] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    # Return
    NNZ <- length(which(PC0[upper.tri(PC0)] != 0))
    cat("- Retained elements: ", NNZ, "\n")
    cat("- Corresponding to", round(NNZ/NR,4) * 100,"% of possible edges \n")
    cat(" \n")

    if (output == "heavy"){
      if (stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        return(PC0)
      }
      if (!stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        colnames(P0)  = rownames(P0)  <- colnames(P)
        return(list(sparseParCor = PC0, sparsePrecision = P0))
      }
    }
    if (output == "light"){
      return(list(zeros = which(PC0 == 0, arr.ind = TRUE),
                  nonzeros = which(PC0 != 0, arr.ind = TRUE)))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Loss/Entropy/Fit Evaluation
##
##------------------------------------------------------------------------------







#' Evaluate regularized precision under various loss functions
#' 
#' Function that evaluates an estimated and possibly regularized precision
#' matrix under various loss functions. The loss functions are formulated in
#' precision terms. This function may be used to estimate the risk (vis-a-vis,
#' say, the true precision matrix) of the various ridge estimators employed.
#' 
#' Let \eqn{\mathbf{\Omega}} denote a generic \eqn{(p \times p)} population
#' precision matrix and let \eqn{\hat{\mathbf{\Omega}}(\lambda)} denote a
#' generic ridge estimator of the precision matrix under generic regularization
#' parameter \eqn{\lambda} (see also \code{\link{ridgeP}}). The function then
#' considers the following loss functions: \enumerate{ \item Squared Frobenius
#' loss, given by: \deqn{ L_{F}[\hat{\mathbf{\Omega}}(\lambda),
#' \mathbf{\Omega}] = \|\hat{\mathbf{\Omega}}(\lambda) -
#' \mathbf{\Omega}\|_{F}^{2}; } \item Quadratic loss, given by: \deqn{
#' L_{Q}[\hat{\mathbf{\Omega}}(\lambda), \mathbf{\Omega}] =
#' \|\hat{\mathbf{\Omega}}(\lambda) \mathbf{\Omega}^{-1} -
#' \mathbf{I}_{p}\|_{F}^{2}.  } } The argument \code{T} is considered to be the
#' true precision matrix when \code{precision = TRUE}. If \code{precision}
#' \code{= FALSE} the argument \code{T} is considered to represent the true
#' covariance matrix. This statement is needed so that the loss is properly
#' evaluated over the precision, i.e., depending on the value of the
#' \code{logical} argument \code{precision} inversions are employed where
#' needed.
#' 
#' The function can be employed to assess the risk of a certain ridge precision
#' estimator (see also \code{\link{ridgeP}}). The risk \eqn{\mathcal{R}_{f}} of
#' the estimator \eqn{\hat{\mathbf{\Omega}}(\lambda)} given a loss function
#' \eqn{L_{f}}, with \eqn{f \in \{F, Q\}} can be defined as the expected loss:
#' \deqn{ \mathcal{R}_{f}[\hat{\mathbf{\Omega}}(\lambda)] =
#' \mathrm{E}\{L_{f}[\hat{\mathbf{\Omega}}(\lambda), \mathbf{\Omega}]\}, }
#' which can be approximated by the mean or median of losses over repeated
#' simulation runs.
#' 
#' @param E Estimated (possibly regularized) precision \code{matrix}.
#' @param T True (population) covariance or precision \code{matrix}.
#' @param precision A \code{logical} indicating if T is a precision matrix.
#' @param type A \code{character} indicating which loss function is to be used.
#' Must be one of: "frobenius", "quadratic".
#' @return Function returns a \code{numeric} representing the loss under the
#' chosen loss function.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}}
#' @references van Wieringen, W.N. & Peeters, C.F.W. (2016).  Ridge Estimation
#' of Inverse Covariance Matrices from High-Dimensional Data, Computational
#' Statistics & Data Analysis, vol. 103: 284-303.  Also available as
#' arXiv:1403.0904v3 [stat.ME].
#' @examples
#' 
#' ## Define population covariance
#' set.seed(333)
#' p = 25
#' n = 1000
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Truecov <- covML(X)
#' 
#' ## Obtain sample
#' samples <- X[sample(nrow(X), 10), ]
#' Cxx <- covML(samples)
#' 
#' ## Obtain regularized precision
#' P <- ridgeP(Cxx, 10, type = "Alt")
#' 
#' ## Evaluate estimated precision against population
#' ## precision under Frobenius loss
#' loss(P, Truecov, precision = FALSE, type = "frobenius")
#' 
#' @export loss
loss <- function(E, T, precision = TRUE, type = c("frobenius", "quadratic")){
  ##############################################################################
  # - Function evualuating various loss functions on the precision
  # - E         > Estimated (possibly regularized) precision matrix
  # - T         > True (population) covariance or precision matrix
  # - precision > Logical indicating if T is a precision matrix (when TRUE)
  # - type      > character indicating which loss function is to be used
  ##############################################################################

  if (!is.matrix(E)){
    stop("Input (E) is of wrong class")
  }
  else if (!isSymmetric(E)){
    stop("E should be a symmetric matrix")
  }
  else if (!is.matrix(T)){
    stop("Input (T) is of wrong class")
  }
  else if (!isSymmetric(T)){
    stop("T should be a symmetric matrix")
  }
  else if (dim(E)[1] != dim(T)[1]){
    stop("E and T should be of the same dimension")
  }
  else if (class(precision) != "logical"){
    stop("Input (precision) is of wrong class")
  }
  else if (missing(type)){
    stop("Need to specify loss type ('frobenius' or 'quadratic')")
  }
  else if (!(type %in% c("frobenius", "quadratic"))){
    stop("type should be one of {'frobenius', 'quadratic'}")
  }
  else {
    # Frobenius loss
    if (type == "frobenius"){
      if (precision)  {loss <- .FrobeniusLoss(E, T)}
      if (!precision) {loss <- .FrobeniusLoss(E, solve(T))}
    }

    # Quadratic loss
    if (type == "quadratic"){
      if (precision)  {loss <- .QuadraticLoss(E, solve(T))}
      if (!precision) {loss <- .QuadraticLoss(E, T)}
    }

    # Return
    return(loss)
  }
}









#' Kullback-Leibler divergence between two multivariate normal distributions
#' 
#' Function calculating the Kullback-Leibler divergence between two
#' multivariate normal distributions.
#' 
#' The Kullback-Leibler (KL) information (Kullback and Leibler, 1951; also
#' known as relative entropy) is a measure of divergence between two
#' probability distributions. Typically, one distribution is taken to represent
#' the `true' distribution and functions as the reference distribution while
#' the other is taken to be an approximation of the true distribution. The
#' criterion then measures the loss of information in approximating the
#' reference distribution. The KL divergence between two \eqn{p}-dimensional
#' multivariate normal distributions
#' \eqn{\mathcal{N}^{0}_{p}(\boldsymbol{\mu}_{0}, \mathbf{\Sigma}_{0})} and
#' \eqn{\mathcal{N}^{1}_{p}(\boldsymbol{\mu}_{1}, \mathbf{\Sigma}_{1})} is
#' given as \deqn{ \mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \| \mathcal{N}^{1}_{p})
#' = \frac{1}{2}\left\{\mathrm{tr}(\mathbf{\Omega}_{1}\mathbf{\Sigma}_{0}) +
#' (\boldsymbol{\mu}_{1} - \boldsymbol{\mu}_{0})^{\mathrm{T}}
#' \mathbf{\Omega}_{1}(\boldsymbol{\mu}_{1} - \boldsymbol{\mu}_{0}) - p -
#' \ln|\mathbf{\Sigma}_{0}| + \ln|\mathbf{\Sigma}_{1}| \right\}, } where
#' \eqn{\mathbf{\Omega} = \mathbf{\Sigma}^{-1}}. The KL divergence is not a
#' proper metric as \eqn{\mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \|
#' \mathcal{N}^{1}_{p}) \neq \mathrm{I}_{KL}(\mathcal{N}^{1}_{p} \|
#' \mathcal{N}^{0}_{p})}. When \code{symmetric = TRUE} the function calculates
#' the symmetric KL divergence (also referred to as Jeffreys information),
#' given as \deqn{ \mathrm{I}_{KL}(\mathcal{N}^{0}_{p} \| \mathcal{N}^{1}_{p})
#' + \mathrm{I}_{KL}(\mathcal{N}^{1}_{p} \| \mathcal{N}^{0}_{p}). }
#' 
#' @param Mtest A \code{numeric} mean vector for the approximating multivariate
#' normal distribution.
#' @param Mref A \code{numeric} mean vector for the true/reference multivariate
#' normal distribution.
#' @param Stest A covariance \code{matrix} for the approximating multivariate
#' normal distribution.
#' @param Sref A covariance \code{matrix} for the true/reference multivariate
#' normal distribution.
#' @param symmetric A \code{logical} indicating if the symmetric version of
#' Kullback-Leibler divergence should be calculated.
#' @return Function returns a \code{numeric} representing the (symmetric)
#' Kullback-Leibler divergence.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}}
#' @references Kullback, S. and Leibler, R.A. (1951). On Information and
#' Sufficiency. Annals of Mathematical Statistics 22: 79-86.
#' @examples
#' 
#' ## Define population
#' set.seed(333)
#' p = 25
#' n = 1000
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cov0  <- covML(X)
#' mean0 <- colMeans(X)
#' 
#' ## Obtain sample from population
#' samples <- X[sample(nrow(X), 10),]
#' Cov1  <- covML(samples)
#' mean1 <- colMeans(samples)
#' 
#' ## Regularize singular Cov1
#' P <- ridgeP(Cov1, 10)
#' CovR <- solve(P)
#' 
#' ## Obtain KL divergence
#' KLdiv(mean1, mean0, CovR, Cov0)
#' 
#' @export KLdiv
KLdiv <- function(Mtest, Mref, Stest, Sref, symmetric = FALSE){
  ##############################################################################
  # - Function that calculates the Kullback-Leibler divergence between two
  #   normal distributions
  # - Mtest     > mean vector approximating m.v. normal distribution
  # - Mref      > mean vector 'true'/reference m.v. normal distribution
  # - Stest     > covariance matrix approximating m.v. normal distribution
  # - Sref      > covariance matrix 'true'/reference m.v. normal distribution
  # - symmetric > logical indicating if original symmetric version of KL div.
  #               should be calculated
  ##############################################################################

  # Dependencies
  # require("base")

  if (class(Mtest) != "numeric"){
    stop("Input (Mtest) is of wrong class")
  }
  else if (class(Mref) != "numeric"){
    stop("Input (Mref) is of wrong class")
  }
  else if (length(Mtest) != length(Mref)){
    stop("Mtest and Mref should be of same length")
  }
  else if (!is.matrix(Stest)){
    stop("Input (Stest) is of wrong class")
  }
  else if (!is.matrix(Sref)){
    stop("Input (Sref) is of wrong class")
  }
  else if (!isSymmetric(Stest)){
    stop("Stest should be symmetric")
  }
  else if (!isSymmetric(Sref)){
    stop("Sref should be symmetric")
  }
  else if (dim(Stest)[1] != length(Mtest)){
    stop("Column and row dimension of Stest should correspond to length Mtest")
  }
  else if (dim(Sref)[1] != length(Mref)){
    stop("Column and row dimension of Sref should correspond to length Mref")
  }
  else if (class(symmetric) != "logical"){
    stop("Input (symmetric) is of wrong class")
  }
  else {
    # Evaluate KL divergence
    KLd <- (sum(diag(solve(Stest) %*% Sref)) +
              t(Mtest - Mref) %*% solve(Stest) %*% (Mtest - Mref) -
              nrow(Sref) - log(det(Sref)) + log(det(Stest)))/2

    # Evaluate (original) symmetric version KL divergence
    if (symmetric){
      KLd <- KLd + (sum(diag(solve(Sref) %*% Stest)) +
                      t(Mref - Mtest) %*% solve(Sref) %*% (Mref - Mtest) -
                      nrow(Sref) - log(det(Stest)) + log(det(Sref)))/2
    }

    # Return
    return(as.numeric(KLd))
  }
}









#' Visual inspection of the fit of a regularized precision matrix
#' 
#' Function aiding the visual inspection of the fit of an estimated (possibly
#' regularized) precision matrix vis-a-vis the sample covariance matrix.
#' 
#' The function outputs various visualizations to aid the visual inspection of
#' an estimated and possibly regularized precision matrix vis-a-vis the sample
#' covariance matrix. The inverse of the estimated precision matrix \code{P} is
#' taken to represent the estimated covariance matrix. The function then
#' outputs a QQ-plot and a heatmap of the observed covariances against the
#' estimated ones. The heatmap has the estimated covariances as
#' lower-triangular elements and the observed covariances as the
#' upper-triangular elements. The function outputs analogous plots for the
#' estimated and observed correlations. In case the observed covariance matrix
#' \code{S} is non-singular also a QQ-plot an a heatmap are generated for the
#' estimated and observed partial correlations.
#' 
#' The function generates files with extension \code{fileType} under default
#' output names. These files are stored in the directory \code{dir} (default is
#' the working directory). To avoid overwriting of files when working in a
#' single directory one may employ the argument \code{nameExt}. By using
#' \code{nameExt} the default output names are extended with a character of
#' choice.
#' 
#' @param Phat (Regularized) estimate of the precision \code{matrix}.
#' @param S Sample covariance \code{matrix}
#' @param diag A \code{logical} determining if the diagonal elements should be
#' retained for plotting.
#' @param fileType A \code{character} determining the output file type. Must be
#' one of: "pdf", "eps".
#' @param nameExt A \code{character} determining the extension of default
#' output names generated by the function.
#' @param dir A \code{character} specifying the directory in which the visual
#' output is stored.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}}
#' @examples
#' 
#' \dontrun{
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain regularized precision matrix
#' P <- ridgeP(Cx, lambda = 10, type = 'Alt')
#' 
#' ## Evaluate visually fit of regularized precision matrix vis-a-vis sample covariance
#' evaluateSfit(P, Cx, diag = FALSE, fileType = "pdf", nameExt = "test")}
#' 
#' @export evaluateSfit
evaluateSfit <- function(Phat, S, diag = FALSE, fileType = "pdf", nameExt = "",
                         dir = getwd()){
  ##############################################################################
  # - Function aiding the visual inspection of the fit of the estimated
  #   (possibly regularized) precision matrix vis-a-vis the sample
  #   covariance matrix
  # - Phat     > (regularized) estimate of the precision matrix
  # - S        > sample covariance matrix
  # - diag     > logical determining treatment diagonal elements for plots
  # - fileType > signifies filetype of output
  # - nameExt  > character giving extension of default output names.
  #              Circumvents overwriting of output when working in single
  #              directory.
  # - dir      > specifies the directory in which the visual output is stored
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")

  if (!is.matrix(Phat)){
    stop("Input (Phat) should be a matrix")
  }
  else if (!isSymmetric(Phat)){
    stop("Input (Phat) should be a symmetric matrix")
  }
  else if (all(diag(Phat) == 1)){
    stop("Input (Phat) should be a nonstandardized precision matrix")
  }
  else if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (all(diag(S) == 1)){
    stop("Input (S) should be a nonstandardized covariance matrix")
  }
  else if (class(diag) != "logical"){
    stop("Input (diag) is of wrong class")
  }
  else if (missing(fileType)){
    stop("Need to specify type of output file ('pdf' or 'eps')")
  }
  else if (!(fileType %in% c("pdf", "eps"))){
    stop("fileType should be one of {'pdf', 'eps'}")
  }
  else if (class(nameExt) != "character"){
    stop("Input (nameExt) is of wrong class")
  }
  else if (class(dir) != "character"){
    stop("Specify directory for output as 'character'")
  }
  else {
    # Obtain estimated covariance matrix
    Shat <- solve(Phat)


    print("Visualizing covariance fit")
    # plot 1: QQ-plot of covariances
    if (fileType == "pdf"){
      pdf(paste(dir, "QQplot_covariances_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "QQplot_covariances_", nameExt, ".eps"))
    }
    if (diag){
      cObs <- as.numeric(S[upper.tri(S, diag = TRUE)])
      cFit <- as.numeric(Shat[upper.tri(Shat, diag = TRUE)])
    }
    if (!diag){
      cObs <- as.numeric(S[upper.tri(S)])
      cFit <- as.numeric(Shat[upper.tri(Shat)])
    }
    op <- par(pty = "s")
    qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample covariances",
           ylab = "fits", main = "QQ-plot, covariances")
    lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          col = "grey", lty = 2)
    par(op); dev.off()

    # plot 2: Comparison of covariances by heatmap
    if (fileType == "pdf"){
      pdf(paste(dir, "heatmap_covariances_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "heatmap_covariances_", nameExt, ".eps"))
    }
    op  <- par(pty = "s")
    slh <- S
    slh[lower.tri(slh)] <- Shat[lower.tri(Shat)]
    gplot <- edgeHeat(slh, diag = diag, legend = FALSE, main = "Covariances")
    print(gplot); par(op); dev.off()


    print("Visualizing correlation fit")
    # plot 3: QQ-plot of correlations
    if (fileType == "pdf"){
      pdf(paste(dir, "QQplot_correlations_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "QQplot_correlations_", nameExt, ".eps"))
    }
    if (diag){
      cObs <- as.numeric(cov2cor(S)[upper.tri(S, diag = TRUE)]);
      cFit <- as.numeric(cov2cor(Shat)[upper.tri(Shat, diag = TRUE)])
    }
    if (!diag){
      cObs <- as.numeric(cov2cor(S)[upper.tri(S)]);
      cFit <- as.numeric(cov2cor(Shat)[upper.tri(Shat)])
    }
    op <- par(pty = "s")
    qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample correlations",
           ylab = "fits", main = "QQ-plot, correlations")
    lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          col = "grey", lty = 2)
    par(op); dev.off()

    # plot 4: Comparison of correlations by heatmap
    if (fileType == "pdf"){
      pdf(paste(dir, "heatmap_correlations_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "heatmap_correlations_", nameExt, ".eps"))
    }
    op  <- par(pty = "s")
    slh <- cov2cor(S)
    slh[lower.tri(slh)] <- cov2cor(Shat)[lower.tri(Shat)]
    gplot <- edgeHeat(slh, diag = diag, legend = FALSE, main = "Correlations")
    print(gplot); par(op); dev.off()


    print("Visualizing partial correlation fit")
    # If the sample covariance matrix non-singular,
    # also evaluate partial correlation fit
    if (evaluateS(S, verbose = FALSE)$posEigen){

      # plot 5: QQ-plot of partial correlations
      if (fileType == "pdf"){
        pdf(paste(dir, "QQplot_partCorrelations_", nameExt, ".pdf"))
      }
      if (fileType == "eps"){
        setEPS();
        postscript(paste(dir, "QQplot_partCorrelations_", nameExt, ".eps"))
      }
      if (diag){
        cObs <- as.numeric(pcor(solve(S))[upper.tri(S)], diag = TRUE);
        cFit <- as.numeric(pcor(Phat)[upper.tri(Phat, diag = TRUE)])}
      if (!diag){
        cObs <- as.numeric(pcor(solve(S))[upper.tri(S)]);
        cFit <- as.numeric(pcor(Phat)[upper.tri(Phat)])
      }
      op <- par(pty = "s")
      qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample partial correlations",
             ylab = "fits", main = "QQ-plot, partial correlations")
      lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
            seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
            col = "grey", lty = 2)
      par(op); dev.off()

      # plot 6: Comparison of partial correlations by heatmap
      if (fileType == "pdf"){
        pdf(paste(dir, "heatmap_partCorrelations_", nameExt, ".pdf"))
      }
      if (fileType == "eps"){
        setEPS();
        postscript(paste(dir, "heatmap_partCorrelations_", nameExt, ".eps"))
      }
      op  <- par(pty = "s")
      slh <- pcor(solve(S))
      slh[lower.tri(slh)] <- pcor(Phat)[lower.tri(Phat)]
      gplot <- edgeHeat(slh, diag = diag, legend = FALSE,
                        main = "Partial correlations")
      print(gplot); par(op); dev.off()

    } else {
      print(paste("sample covariance matrix is singular:",
                  "partial correlation fit not visualized"))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Visualization
##
##------------------------------------------------------------------------------







#' Visualize the regularization path
#' 
#' Function that visualizes the regularization paths of the nonredundant
#' elements of a regularized precision matrix against the (range of the)
#' penalty parameter.
#' 
#' The function visualizes the regularization path of the individual elements
#' of a regularized precision matrix against the penalty parameter. The range
#' of the penalty parameter is given by [\code{lambdaMin},\code{lambdaMax}].
#' The penalty parameter must be positive such that \code{lambdaMin} must be a
#' positive scalar. The maximum allowable value of \code{lambdaMax} depends on
#' the type of ridge estimator employed. For details on the type of ridge
#' estimator one may use (one of: "Alt", "ArchI", "ArchII") see
#' \code{\link{ridgeP}}.
#' 
#' Regularization paths may be visualized for (partial) correlations,
#' covariances and precision elements. The type of element for which a
#' visualization of the regularization paths is desired can be indicated by the
#' argument \code{plotType}. When \code{vertical = TRUE} a vertical line is
#' added at the constant \code{value}. This option can be used to assess
#' whereabouts the optimal penalty obtained by, e.g., the routines
#' \code{\link{optPenalty.LOOCV}} or \code{\link{optPenalty.aLOOCV}}, finds
#' itself along the regularization path.
#' 
#' @param S Sample covariance \code{matrix}.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param step An \code{integer} determining the number of steps in moving
#' through the grid [\code{lambdaMin}, \code{lambdaMax}].
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param plotType A \code{character} indicating the type of element for which
#' a visualization of the regularization paths is desired. Must be one of:
#' "pcor", "cor", "cov", "prec".
#' @param diag A \code{logical} indicating if the diagonal elements should be
#' retained for visualization.
#' @param vertical A \code{logical} indicating if output graph should come with
#' a vertical line at a pre-specified value for the penalty parameter.
#' @param value A \code{numeric} indicating a pre-specified value for the
#' penalty parameter.
#' @param verbose A \code{logical} indicating if information on progress should
#' be printed on screen.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}},
#' \code{\link{optPenalty.LOOCV}}, \code{\link{optPenalty.aLOOCV}},
#' \code{\link{default.target}}
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Visualize regularization paths
#' ridgePathS(Cx, .001, 50, 200, plotType = "pcor")
#' 
#' @export ridgePathS
ridgePathS <- function (S, lambdaMin, lambdaMax, step, type = "Alt",
                        target = default.target(S), plotType = "pcor",
                        diag = FALSE, vertical = FALSE, value, verbose = TRUE){
  ##############################################################################
  # - Function that visualizes the regularization path
  # - Regularization path may be visualized for (partial) correlations,
  #   covariances and precision elements
  # - S         > sample covariance/correlation matrix
  # - lambdaMin > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax > maximum value penalty parameter (dependent on 'type')
  # - step      > determines the coarseness in searching the grid
  #               [lambdaMin, lambdaMax]
  # - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target    > target (precision terms) for Type I estimators,
  #               default = default.target(S)
  # - plotType  > specificies the elements for which the regularization path is
  #               to be visualized.
  #               Must be one of {"pcor", "cor", "cov", "prec"},
  #               default = "pcor"
  # - diag      > logical indicating if the diagonal elements should be retained
  #               for plotting, default = FALSE.
  # - vertical  > optional argument for visualization vertical line in graph
  #               output, default = FALSE
  #               Can be used to indicate the value of, e.g., the optimal
  #               penalty as indicated by some
  #               routine. Can be used to assess the whereabouts of this optimal
  #               penalty along the regularization path.
  # - value     > indicates constant on which to base vertical line when
  #               vertical = TRUE
  # - verbose   > logical indicating if intermediate output should be printed
  #               on screen
  ##############################################################################
  # Dependencies
  # require("base")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(S)){
    stop("input (S) should be a matrix")
  }
  if (!isSymmetric(S)){
    stop("Input (S) should be a covariance matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("lambdaMax must be larger than lambdaMin")
  }
  else if (class(step) != "numeric") {
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (class(plotType) != "character") {
    stop("Input (plotType) is of wrong class")
  }
  else if (!(plotType %in% c("pcor", "cor", "cov", "prec"))){
    stop("Input (plotType) should be one of {'pcor', 'cor', 'cov', 'prec'}")
  }
  else if (length(nchar(plotType)) != 1){
    stop("Input (plotType) should be exactly one of {'pcor', 'cor', 'cov', ",
         "'prec'}")
  }
  if (class(diag) != "logical") {
    stop("Input (diag) is of wrong class")
  }
  else if (class(vertical) != "logical"){
    stop("Input (vertical) is of wrong class")
  }
  else {
    # Set preliminaries
    lambdas  <- seq(lambdaMin, lambdaMax, len = step)
    YforPlot <- numeric()

    # Calculate paths
    if (verbose){cat("Calculating...", "\n")}
    if (type == "Alt" & all(target == 0)){
      if (!isSymmetric(target)){
        stop("Input (target) should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Inputs ('S' and 'target') should be of the same dimension")
      } else {
        Spectral <- eigen(S, symmetric = TRUE)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k])
          P         <- Spectral$vectors %*% diag(1/Eigshrink) %*% t(Spectral$vectors)
          if (plotType=="pcor"){
            YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
          }
          if (plotType=="prec"){
            YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
          }
          if (plotType=="cov") {
            YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
          }
          if (plotType=="cor") {
            YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
          }
          if (verbose){
            cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
          }
        }
      }
    } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
      if (!isSymmetric(target)){
        stop("Input (target) should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Inputs ('S' and 'target') should be of the same dimension")
      } else if (any(diag(target) <= 0)){
        stop("Input (target) should be p.d.")
      } else {
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k], const = varPhi)
          P         <- Spectral$vectors %*% diag(1/Eigshrink) %*% t(Spectral$vectors)
          if (plotType=="pcor"){
            YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
          }
          if (plotType=="prec"){
            YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
          }
          if (plotType=="cov") {
            YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
          }
          if (plotType=="cor") {
            YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
          }
          if (verbose){
            cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
          }
        }
      }
    } else {
      for (k in 1:length(lambdas)){
        P <- ridgeP(S, lambdas[k], type = type, target = target)
        if (plotType=="pcor"){
          YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
        }
        if (plotType=="prec"){
          YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
        }
        if (plotType=="cov") {
          YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
        }
        if (plotType=="cor") {
          YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
        }
        if (verbose){
          cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
        }
      }
    }

    # Visualize
    if (plotType=="cor") {ylabel <- "penalized correlation"}
    if (plotType=="cov") {ylabel <- "penalized covariances"}
    if (plotType=="pcor"){ylabel <- "penalized partial correlation"}
    if (plotType=="prec"){ylabel <- "penalized precision elements"}
    if (type == "Alt"){Main = "Alternative ridge estimator"}
    if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
    if (type == "ArchII"){Main = "Archetypal II ridge estimator"}

    plot(YforPlot[1,] ~ log(lambdas), axes = FALSE, xlab = "ln(penalty value)",
         ylab = ylabel, main = Main, col = "white",
         ylim = c(min(YforPlot), max(YforPlot)))
    for (k in 1:nrow(YforPlot)){
      lines(YforPlot[k, ] ~ log(lambdas), col = k, lty = k)
    }
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    if (vertical){
      if (missing(value)){
        stop("Need to specify input (value)")
      } else if (class(value) != "numeric"){
        stop("Input (value) is of wrong class")
      } else if (length(value) != 1){
        stop("Input (value) must be a scalar")
      } else if (value <= 0){
        stop("Input (value) must be positive")
      } else {
        abline(v = log(value), col = "red")
      }
    }
  }
}



if (getRversion() >= "2.15.1") utils::globalVariables(c("X1", "X2", "value"))







#' Visualize (precision) matrix as a heatmap
#' 
#' Function that visualizes a (precision) matrix as a heatmap. May be used to
#' assess visually the elements of a single (possibly sparsified precision)
#' matrix. May also be used in assessing the performance of edge selection
#' techniques.
#' 
#' This function utilizes
#' \href{https://cran.r-project.org/package=ggplot2ggplot2} (Wickham, 2009) to
#' visualize a matrix as a heatmap: a false color plot in which the individual
#' matrix entries are represented by colors. \code{lowColor} determines the
#' color scale for matrix entries in the negative range. \code{highColor}
#' determines the color scale for matrix entries in the positive range. For the
#' colors supported by the arguments \code{lowColor} and \code{highColor}, see
#' \url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}. White entries in
#' the plot represent the midscale value of 0. One can opt to set the diagonal
#' entries to the midscale color of white when one is interested in
#' (heatmapping) the off-diagonal elements only. To achieve this, set
#' \code{diag = FALSE}. Naturally, the \code{diag} argument is only used when
#' the input matrix \code{M} is a square matrix.
#' 
#' The intended use of the function is to visualize a, possibly sparsified,
#' precision matrix as a heatmap. The function may also be used, in a graphical
#' modeling setting, to assess the performance of edge selection techniques.
#' However, the function is quite general, in the sense that it can represent
#' any \code{matrix} as a heatmap.
#' 
#' @param M (Possibly sparsified precision) \code{matrix}.
#' @param lowColor A \code{character} that determines the color scale in the
#' negative range.
#' @param highColor A \code{character} that determines the color scale in the
#' positive range.
#' @param textsize A \code{numeric} scaling the text size of row and column
#' labels.
#' @param diag A \code{logical} determining if the diagonal elements of the
#' matrix should be included in the color scaling. This argument is only used
#' when \code{M} is a square \code{matrix}.
#' @param legend A \code{logical} indicating whether a color legend should be
#' included.
#' @param main A \code{character} giving the main figure title.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{covML}}, \code{\link{ridgeP}}, \code{\link{sparsify}}
#' @references Wickham, H. (2009). ggplot2: elegant graphics for data analysis.
#' New York: Springer.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain regularized precision matrix
#' P <- ridgeP(Cx, lambda = 10, type = "Alt")
#' 
#' ## Obtain sparsified partial correlation matrix
#' PC0 <- sparsify(P, threshold = "localFDR", FDRcut = .8)$sparseParCor
#' 
#' ## Visualize sparsified partial correlation matrix as heatmap
#' edgeHeat(PC0)
#' 
#' @export edgeHeat
edgeHeat <- function(M, lowColor = "blue", highColor = "red", textsize = 10,
                     diag = TRUE, legend = TRUE, main = ""){
  ##############################################################################
  # - function that visualizes precision matrix as a heatmap
  # - can be used to assess (visually) the performance of set of graphical
  #   modeling techniques
  # - M         > Precision matrix
  # - lowColor  > determines color scale in the negative range, default = "blue"
  # - highColor > determines color scale in the positive range, default = "red"
  # - textsize  > set textsize row and column labels, default = 10
  # - diag      > logical determining treatment diagonal elements M. If FALSE,
  #               then the diagonal elements are given the midscale color of
  #               white; only when M is a square matrix
  # - legend    > optional inclusion of color legend, default = TRUE
  # - main      > character specifying the main title, default = ""
  ##############################################################################

  # Dependencies
  #require("ggplot2")
  #require("reshape")

  if (!is.matrix(M)){
    stop("Supply 'M' as matrix")
  }
  else if (class(lowColor) != "character"){
    stop("Input (lowColor) is of wrong class")
  }
  else if (length(lowColor) != 1){
    stop("Length lowColor must be one")
  }
  else if (class(highColor) != "character"){
    stop("Input (highColor) is of wrong class")
  }
  else if (length(highColor) != 1){
    stop("Length highColor must be one")
  }
  else if (class(textsize) != "numeric"){
    stop("Input (textsize) is of wrong class")
  }
  else if (length(textsize) != 1){
    stop("Length textsize must be one")
  }
  else if (textsize <= 0){
    stop("textsize must be positive")
  }
  else if (class(diag) != "logical"){
    stop("Input (diag) is of wrong class")
  }
  else if (class(legend) != "logical"){
    stop("Input (legend) is of wrong class")
  }
  else if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  else {
    # Put matrix in data format
    if (nrow(M) == ncol(M) & !diag) {diag(M) <- 0}
    Mmelt    <- melt(M)
    Mmelt$X1 <-
      factor(as.character(Mmelt$X1), levels = unique(Mmelt$X1), ordered = TRUE)
    Mmelt$X2 <-
      factor(as.character(Mmelt$X2), levels = unique(Mmelt$X2), ordered = TRUE)

    # Visualize
    if (legend){
      ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() +
        scale_fill_gradient2("", low = lowColor,  mid = "white",
                             high = highColor, midpoint = 0) +
        theme(axis.ticks = element_blank()) +
        theme(axis.text.y = element_text(size = textsize)) +
        theme(axis.text.x = element_text(angle = -90, vjust = .5,
                                         hjust = 0, size = textsize)) +
        xlab(" ") + ylab(" ") +
        ylim(rev(levels(Mmelt$X1))) +
        ggtitle(main)
    } else {
      ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() +
        scale_fill_gradient2("", low = lowColor,  mid = "white",
                             high = highColor, midpoint = 0) +
        theme(axis.ticks = element_blank()) +
        theme(axis.text.y = element_text(size = textsize)) +
        theme(axis.text.x = element_text(angle = -90, vjust = .5,
                                         hjust = 0, size = textsize)) +
        xlab(" ") + ylab(" ") +
        ylim(rev(levels(Mmelt$X1))) +
        ggtitle(main) +
        theme(legend.position = "none")
    }
  }
}









#' Visualize undirected graph
#' 
#' Function that visualizes the sparsified precision matrix as an undirected
#' graph.
#' 
#' The intended use of this function is to visualize a sparsified
#' precision/partial correlation matrix as an undirected graph. When \code{type
#' = "plain"} a plain undirected graph is given representing the conditional
#' (in)dependencies exemplified by the sparsified precision.
#' 
#' When \code{type = "fancy"} a more elaborate graph is given in which dashed
#' lines indicate negative partial correlations while solid lines indicate
#' positive partial correlations, and in which grey lines indicate strong
#' edges. Strong edges are deemed such by setting \code{cut}. If a the absolute
#' value of a precision element \eqn{\geq} \code{cut} the corresponding edge is
#' deemed strong and colored grey in the graph. The argument \code{cut} is thus
#' only used when \code{type = "fancy"}.
#' 
#' When \code{type = "weighted"} an undirected graph is given in which edge
#' thickness represents the strength of the partial correlations. The
#' \code{nEcolor} colored edges then represent negative partial correlations
#' while \code{pEcolor} colored edges represent positive partial correlations.
#' (Relative) edge thickness in this type of graph can be set by the argument
#' \code{scale}. The arguments \code{scale}, \code{nEcolor}, and \code{pEcolor}
#' are thus only used when \code{type = "weighted"}.
#' 
#' The default layout gives a circular placement of the vertices. Most layout
#' functions supported by \code{\link{igraph}} are supported (the function is
#' partly a wrapper around certain \code{\link{igraph}} functions). The igraph
#' layouts can be invoked by a \code{character} that mimicks a call to a
#' \code{\link{igraph}} layout functions in the \code{lay} argument. When using
#' \code{lay = NULL} one can specify the placement of vertices with the
#' \code{coords} argument. The row dimension of this matrix should equal the
#' number of (pruned) vertices. The column dimension then should equal 2 (for
#' 2D layouts) or 3 (for 3D layouts). The \code{coords} argument can also be
#' viewed as a convenience argument as it enables one, e.g., to layout a graph
#' according to the coordinates of a previous call to \code{Ugraph}. If both
#' the the lay and the coords arguments are not \code{NULL}, the lay argument
#' takes precedence
#' 
#' The legend allows one to specify the kind of variable the vertices
#' represent, such as, e.g., mRNA transcripts. The arguments \code{label},
#' \code{Lcex}, and \code{PTcex} are only used when \code{legend = TRUE}.
#' 
#' If \code{prune = TRUE} the vertices of degree 0 (vertices not implicated by
#' any edge) are removed. For the colors supported by the arguments
#' \code{Vcolor}, \code{VBcolor}, \code{VLcolor}, \code{pEcolor}, and
#' \code{nEcolor} see \url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}.
#' 
#' @param M (Possibly sparsified) precision \code{matrix}
#' @param type A \code{character} indicating the type of graph to be produced.
#' Must be one of: "plain", "fancy", "weighted".
#' @param lay A \code{character} mimicking a call to \code{\link{igraph}}
#' layout functions. Determines the placement of vertices.
#' @param coords A \code{matrix} containing coordinates. Alternative to the
#' lay-argument for determining the placement of vertices.
#' @param Vsize A \code{numeric} determining the vertex size.
#' @param Vcex A \code{numeric} determining the size of the vertex labels.
#' @param Vcolor A \code{character} (scalar or vector) determining the vertex
#' color.
#' @param VBcolor A \code{character} determining the color of the vertex
#' border.
#' @param VLcolor A \code{character} determining the color of the vertex
#' labels.
#' @param prune A \code{logical} determining if vertices of degree 0 should be
#' removed.
#' @param legend A \code{logical} indicating if the graph should come with a
#' legend.
#' @param label A \code{character} giving a name to the legend label.
#' @param Lcex A \code{numeric} determining the size of the legend box.
#' @param PTcex A \code{numeric} determining the size of the exemplary vertex
#' in the legend box.
#' @param cut A \code{numeric} indicating the cut-off for indicating strong
#' edges when \code{type = "fancy"}.
#' @param scale A \code{numeric} representing a scale factor for visualizing
#' strength of edges when \code{type = "weighted"}.
#' @param pEcolor A \code{character} determining the color of the edges tied to
#' positive precision elements. Only when \code{type = "weighted"}.
#' @param nEcolor A \code{character} determining the color of the edges tied to
#' negative precision elements. Only when \code{type = "weighted"}.
#' @param main A \code{character} giving the main figure title.
#' @return The function returns a graph. The function also returns a
#' \code{matrix} object containing the coordinates of the vertices in the given
#' graph.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.LOOCV}},
#' \code{\link{optPenalty.aLOOCV}}, \code{\link{sparsify}}
#' @references Csardi, G. and Nepusz, T. (2006). The igraph software package
#' for complex network research. InterJournal, Complex Systems 1695.
#' http://igraph.sf.net
#' 
#' van Wieringen, W.N. & Peeters, C.F.W. (2016). Ridge Estimation of Inverse
#' Covariance Matrices from High-Dimensional Data, Computational Statistics &
#' Data Analysis, vol. 103: 284-303. Also available as arXiv:1403.0904v3
#' [stat.ME].
#' 
#' van Wieringen, W.N. & Peeters, C.F.W. (2015). Application of a New Ridge
#' Estimator of the Inverse Covariance Matrix to the Reconstruction of
#' Gene-Gene Interaction Networks. In: di Serio, C., Lio, P., Nonis, A., and
#' Tagliaferri, R. (Eds.) `Computational Intelligence Methods for
#' Bioinformatics and Biostatistics'. Lecture Notes in Computer Science, vol.
#' 8623. Springer, pp. 170-179.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)
#' 
#' ## Determine support regularized standardized precision under optimal penalty
#' PC0 <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor
#' 
#' ## Obtain graphical representation
#' Ugraph(PC0, type = "fancy", cut = 0.07)
#' 
#' ## Obtain graphical representation with Fruchterman-Reingold layout
#' Ugraph(PC0, type = "fancy", lay = "layout_with_fr", cut = 0.07)
#' 
#' ## Add pruning
#' Ugraph(PC0, type = "fancy", lay = "layout_with_fr",
#'        cut = 0.07, prune = TRUE)
#' 
#' ## Obtain graph and its coordinates
#' Coordinates <- Ugraph(PC0, type = "fancy", lay = "layout_with_fr",
#'                       cut = 0.07, prune = TRUE)
#' Coordinates
#' 
#' @export Ugraph
Ugraph <- function(M, type = c("plain", "fancy", "weighted"),
                   lay = "layout_in_circle", coords = NULL, Vsize = 15,
                   Vcex = 1, Vcolor = "orangered", VBcolor = "darkred",
                   VLcolor = "black", prune = FALSE, legend = FALSE,
                   label = "", Lcex = 1.3, PTcex = 4, cut = .5,
                   scale = 10, pEcolor = "black", nEcolor = "grey",
                   main = ""){
  ##############################################################################
  # - Function that visualizes the sparsified precision matrix as an undirected
  #   graph
  # - Function is partly a wrapper around certain 'igraph' functions
  # - M       > (Possibly sparsified) precision matrix
  # - type    > graph type: 'plain' gives plain undirected graph. 'fancy' gives
  #             undirected graph in which dashed lines indicate negative partial
  #             correlations while solid lines indicate positive partial
  #             correlations, and in which grey lines indicate strong edges.
  #             'weighted' gives an undirected graph in which edge thickness
  #             indicates the strenght of the partial correlations. Grey lines
  #             then indicate negative partial correlations while black lines
  #             represent positive partial correlations.
  # - lay     > determines layout of the graph. Most layouts in 'layout{igraph}'
  #             are accepted. Default = layout_in_circle.
  # - coords  > matrix of coordinates to determine layout of the graph.
  #             The row dimension should equal the number of (pruned) vertices.
  #             The column dimension should equal 2 (for 2D layouts) or
  #             3 (for 3D layouts). Enables one, e.g., to layout the graph
  #             according to the coordinates of a previous call to Ugraph.
  #             If both the the lay and the coords arguments are not NULL,
  #             the lay argument takes precedence
  # - Vsize   > gives vertex size, default = 15
  # - Vcex    > gives size vertex labels, default = 1
  # - Vcolor  > gives vertex color, default = "orangered", must be character.
  #             May also be a character vector
  # - VBcolor > gives color of the vertex border, default = "darkred"
  # - VLcolor > gives color of the vertex labels, default = "black"
  # - prune   > logical indicating if vertices of degree 0 should be removed
  # - legend  > optional inclusion of color legend, default = FALSE
  # - label   > character label for the endogenous variables, default = "";
  #             only when legend = TRUE
  # - Lcex    > scaling legend box, default = 1.3; only when legend = TRUE
  # - PTcex   > scaling node in legend box, default = 4; only when legend = TRUE
  # - cut     > cut-off for indication of strong edge, default = .5; only when
  #             type = "fancy"
  # - scale   > scale factor for visualizing strenght of edges, default = 10;
  #             only when type = "weighted"
  # - pEcolor > gives edge color for edges tied to positive precision elements,
  #             default = "black"; only when type = "weighted"
  # - nEcolor > gives edge color for edges tied to negative precision elements,
  #             default = "grey"; only when type = "weighted"
  # - main    > character specifying heading figure, default = ""
  ##############################################################################

  # Dependencies
  # require("igraph")
  # require("reshape")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be square matrix")
  }
  else if (missing(type)){
    stop("Need to specify graph type ('plain' or 'fancy' or 'weighted')")
  }
  else if (!(type %in% c("plain", "fancy", "weighted"))){
    stop("type should be one of {'plain', 'fancy', 'weighted'}")
  }
  else if (!((length(intersect(lay,c("layout_as_star", "layout_as_tree",
                                     "layout_in_circle", "layout_nicely",
                                     "layout_with_dh", "layout_with_gem",
                                     "layout_with_graphopt", "layout_on_grid",
                                     "layout_with_mds", "layout_components",
                                     "layout_on_sphere", "layout_randomly",
                                     "layout_with_fr", "layout_with_kk",
                                     "layout_with_lgl"))) > 0) | is.null(lay))){
    stop("lay should be 'NULL' or one of
         {'layout_as_star', 'layout_as_tree',
         'layout_in_circle', 'layout_nicely',
         'layout_with_dh', 'layout_with_gem',
         'layout_with_graphopt', 'layout_on_grid',
         'layout_with_mds', 'layout_components',
         'layout_on_sphere', 'layout_randomly',
         'layout_with_fr', 'layout_with_kk',
         'layout_with_lgl'}")
  }
  else if (!is.null(coords) & !inherits(coords, "matrix")){
    stop("Input (coords) is of wrong class")
  }
  else if (is.null(lay) & is.null(coords)){
    stop("Input (lay) and input (coords) cannot be both NULL")
  }
  else if (class(Vsize) != "numeric"){
    stop("Input (Vsize) is of wrong class")
  }
  else if (length(Vsize) != 1){
    stop("Length Vsize must be one")
  }
  else if (Vsize <= 0){
    stop("Vsize must be positive")
  }
  else if (class(Vcex) != "numeric"){
    stop("Input (Vcex) is of wrong class")
  }
  else if (length(Vcex) != 1){
    stop("Length Vcex must be one")
  }
  else if (Vcex <= 0){
    stop("Vcex must be positive")
  }
  else if (class(Vcolor) != "character"){
    stop("Input (Vcolor) is of wrong class")
  }
  else if (length(Vcolor) != 1 & length(Vcolor) != nrow(M)){
    stop("Length Vcolor must be either one
         or equal to row (or column) dimension of M")
  }
  else if (class(VBcolor) != "character"){
    stop("Input (VBcolor) is of wrong class")
  }
  else if (length(VBcolor) != 1){
    stop("Length VBcolor must be one")
  }
  else if (class(VLcolor) != "character"){
    stop("Input (VLcolor) is of wrong class")
  }
  else if (length(VLcolor) != 1){
    stop("Length VLcolor must be one")
  }
  else if (class(prune) != "logical"){
    stop("Input (prune) is of wrong class")
  }
  else if (class(legend) != "logical"){
    stop("Input (legend) is of wrong class")
  }
  else if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  else {
    # Preliminaries
    AM <- adjacentMat(M)
    GA <- graph.adjacency(AM, mode = "undirected")
    if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}

    # Layout specification
    if(is.null(lay)){
      if(dim(coords)[1] != length(V(GA))){
        stop("Row dimension of input (coords) does not match the
             number of vertices to be plotted")
      } else if (dim(coords)[2] > 3){
        stop("Column dimension of input (coords) exceeds the number
             of dimensions that can be visualized")
      } else {lays = coords}
      }
    else{
      if(lay == "layout_as_star"){
        lays = igraph::layout_as_star(GA)}
      if(lay == "layout_as_tree")
      {lays = igraph::layout_as_tree(GA)}
      if(lay == "layout_in_circle"){
        lays = igraph::layout_in_circle(GA)}
      if(lay == "layout_nicely"){
        lays = igraph::layout_nicely(GA)}
      if(lay == "layout_with_dh"){
        lays = igraph::layout_with_dh(GA)}
      if(lay == "layout_with_gem"){
        lays = igraph::layout_with_gem(GA)}
      if(lay == "layout_with_graphopt"){
        lays = igraph::layout_with_graphopt(GA)}
      if(lay == "layout_on_grid"){
        lays = igraph::layout_on_grid(GA)}
      if(lay == "layout_with_mds"){
        lays = igraph::layout_with_mds(GA)}
      if(lay == "layout_components"){
        lays = igraph::layout_components(GA)}
      if(lay == "layout_on_sphere"){
        lays = igraph::layout_on_sphere(GA)}
      if(lay == "layout_randomly"){
        lays = igraph::layout_randomly(GA)}
      if(lay == "layout_with_fr"){
        lays = igraph::layout_with_fr(GA)}
      if(lay == "layout_with_kk"){
        lays = igraph::layout_with_kk(GA)}
      if(lay == "layout_with_lgl"){
        lays = igraph::layout_with_lgl(GA)}
    }

    # Plain graph
    if (type == "plain"){
      plot(GA, layout = lays, vertex.size = Vsize, vertex.label.family = "sans",
           vertex.label.cex = Vcex, vertex.color = Vcolor,
           vertex.frame.color = VBcolor,
           vertex.label.color = VLcolor, main = main)
    }

    # Fancy graph
    if (type == "fancy"){
      if (class(cut) != "numeric"){
        stop("Input (cut) is of wrong class")
      } else if (length(cut) != 1){
        stop("Length cut must be one")
      } else if (cut <= 0){
        stop("cut must be positive")
      } else {
        Names <- colnames(M)
        colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
        Mmelt <- melt(M)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]
        E(GA)$weight <- Mmelt$value
        E(GA)$color  <- "black"
        E(GA)[E(GA)$weight < 0]$style <- "dashed"
        E(GA)[E(GA)$weight > 0]$style <- "solid"
        E(GA)[abs(E(GA)$weight) > cut]$color <- "grey"
        plot(GA, layout = lays, vertex.size = Vsize,
             vertex.label.family = "sans", vertex.label.cex = Vcex,
             vertex.color = Vcolor, vertex.frame.color = VBcolor,
             vertex.label.color = VLcolor,
             edge.color = E(GA)$color, edge.lty = E(GA)$style, main = main)
      }
    }

    # Weighted graph
    if (type == "weighted"){
      if (class(scale) != "numeric"){
        stop("Input (scale) is of wrong class")
      } else if (length(scale) != 1){
        stop("Length scale must be one")
      } else if (scale <= 0){
        stop("scale must be positive")
      } else if (class(pEcolor) != "character"){
        stop("Input (pEcolor) is of wrong class")
      } else if (length(pEcolor) != 1){
        stop("Length pEcolor must be one")
      } else if (class(nEcolor) != "character"){
        stop("Input (nEcolor) is of wrong class")
      } else if (length(nEcolor) != 1){
        stop("Length nEcolor must be one")
      } else {
        Names <- colnames(M)
        colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
        Mmelt <- melt(M)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]
        E(GA)$weight <- Mmelt$value
        E(GA)[E(GA)$weight < 0]$color <- nEcolor
        E(GA)[E(GA)$weight > 0]$color <- pEcolor
        plot(GA, layout = lays, vertex.size = Vsize,
             vertex.label.family = "sans", vertex.label.cex = Vcex,
             vertex.color = Vcolor, vertex.frame.color = VBcolor,
             vertex.label.color = VLcolor,
             edge.color = E(GA)$color, edge.width = scale*abs(E(GA)$weight),
             main = main)
      }
    }

    # Legend
    if (legend){
      if (class(label) != "character"){
        stop("Input (label) is of wrong class")
      } else if (length(label) != 1){
        stop("Length label must be one")
      } else if (class(Lcex) != "numeric"){
        stop("Input (Lcex) is of wrong class")
      } else if (length(Lcex) != 1){
        stop("Length Lcex must be one")
      } else if (Lcex <= 0){
        stop("Lcex must be positive")
      } else if (class(PTcex) != "numeric"){
        stop("Input (PTcex) is of wrong class")
      } else if (length(PTcex) != 1){
        stop("Length PTcex must be one")
      } else if (PTcex <= 0){
        stop("PTcex must be positive")
      } else{
        legend("bottomright", label, pch = 20, col = Vcolor,
               cex = Lcex, pt.cex = PTcex)
      }
    }

    # Return
    return(coordinates <- lays)
    }
  }




##------------------------------------------------------------------------------
##
## Functions for Topology Statistics
##
##------------------------------------------------------------------------------







#' Gaussian graphical model network statistics
#' 
#' Function that calculates various network statistics from a sparse precision
#' matrix. The sparse precision matrix is taken to represent the conditional
#' indepence graph of a Gaussian graphical model.
#' 
#' The function calculates various network statistics from a sparse matrix. The
#' input matrix \code{P} is assumed to be a sparse precision or partial
#' correlation matrix. The sparse matrix is taken to represent a conditional
#' independence graph. In the Gaussian setting, conditional independence
#' corresponds to zero entries in the (standardized) precision matrix. Each
#' node in the graph represents a Gaussian variable, and each undirected edge
#' represents conditional dependence in the sense of a nonzero corresponding
#' precision entry.
#' 
#' The function calculates various measures of centrality: node degree,
#' betweenness centrality, closeness centrality, and eigenvalue centrality. It
#' also calculates the number of positive and the number of negative edges for
#' each node. In addition, for each variate the mutual information (with all
#' other variates), the variance, and the partial variance is represented. It
#' is also indicated if the graph is chordal (i.e., triangulated). For more
#' information on network measures, consult, e.g., Newman (2010).
#' 
#' @param sparseP Sparse precision/partial correlation \code{matrix}.
#' @param as.table A \code{logical} indicating if the output should be in
#' tabular format.
#' @return An object of class \code{list} when \code{as.table = FALSE}:
#' \item{degree}{A \code{numeric} vector with the node degree for each node.}
#' \item{betweenness}{A \code{numeric} vector representing the betweenness
#' centrality for each node.} \item{closeness}{A \code{numeric} vector
#' representing the closeness centrality for each node.}
#' \item{eigenCentrality}{A \code{numeric} vector representing the eigenvalue
#' centrality for each node.} \item{nNeg}{An \code{integer} vector representing
#' the number of negative edges for each node.} \item{nPos}{An \code{integer}
#' vector representing the number of positive edges for each node.}
#' \item{chordal}{A \code{logical} indicating if the implied graph is chordal.}
#' \item{mutualInfo}{A \code{numeric} vector with the mutual information (with
#' all other nodes) for each node.} \item{variance}{A \code{numeric} vector
#' representing the variance of each node.} \item{partialVariance}{A
#' \code{numeric} vector representing the partial variance of each node.} When
#' \code{as.table = TRUE} the list items above (with the exception of
#' \code{chordal}) are represented in tabular form as an object of class
#' \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{covML}}, \code{\link{sparsify}},
#' \code{\link{Ugraph}}
#' @references Newman, M.E.J. (2010). "Networks: an introduction", Oxford
#' University Press.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' Cx <- covML(X)
#' 
#' ## Obtain sparsified partial correlation matrix
#' Pridge   <- ridgeP(Cx, 10, type = "Alt")
#' PCsparse <- sparsify(Pridge , threshold = "top")$sparseParCor
#' 
#' ## Represent the graph and calculate GGM network statistics
#' Ugraph(PCsparse, "fancy")
#' \dontrun{GGMnetworkStats(PCsparse)}
#' 
#' @export GGMnetworkStats
GGMnetworkStats <- function(sparseP, as.table = FALSE){
  ##############################################################################
  # - Function that calculates various network statistics from a sparse matrix
  # - Input matrix is assumed to be a sparse precision of partial correlation
  #   matrix
  # - The sparse precision matrix is taken to represent a conditional
  #   independence graph
  # - sparseP  > sparse precision/partial correlation matrix
  # - as.table > logical indicating if output should be returned as table;
  #              default = FALSE
  #
  # - NOTES (network statistics produced):
  # - Node degree
  # - Betweenness centrality
  # - Closeness centrality
  # - Eigenvalue centrality
  # - Number of negative edges for each node
  # - Number of positive edges for each node
  # - Assessment if network/graph is chordal (triangulated)
  # - Mutual information of each variate with all other variates
  # - Variance of each variate (based on inverse sparsified precision matrix)
  # - Partial variance of each variate (= 1 when input matrix is partial
  #   correlation matrix)
  # - Future versions of this function may include additional statistics
  #
  # - REFERENCE:
  # - Newman, M.E.J. (2010), "Networks: an introduction",
  #   Oxford University Press
  ##############################################################################

  # Dependencies
  # require("base")
  # require("igraph")

  if (!is.matrix(sparseP)){
    stop("Input (sparseP) should be a matrix")
  }
  else if (!isSymmetric(sparseP)){
    stop("Input (sparseP) should be a symmetric matrix")
  }
  else if (class(as.table) != "logical"){
    stop("Input (as.table) is of wrong class")
  }
  else{
    # Some warnings
    if (all(sparseP != 0)){
      warning("Given input (sparseP) implies a saturated conditional ",
              "independence graph")
    }
    if (all(sparseP[!diag(nrow(sparseP))] == 0)){
      warning("Given input (sparseP) implies an empty conditional ",
              "independence graph")
    }

    # Obtain corresponding sample covariance matrix
    pvars <- 1/diag(sparseP)
    S     <- solve(sparseP)

    # Calculate nodes' mutual information
    MI <-
      unlist(lapply(1:nrow(S),
                    function(j, S){
                      log(det(S[-j,-j])) -
                        log(det(S[-j,-j] - S[-j,j,drop=FALSE] %*%
                                  S[j,-j,drop=FALSE]/S[j,j]))
                    }, S = S))
    names(MI) <- colnames(sparseP)

    # Signs of edges
    diag(sparseP) <- 0
    nPos <- apply(sign(sparseP), 2, function(Z){ sum(Z == 1) })
    nNeg <- apply(sign(sparseP), 2, function(Z){ sum(Z == -1) })

    # Adjacency to graphical object
    AM  <- adjacentMat(sparseP)
    CIG <- graph.adjacency(AM, mode = "undirected")

    # Return
    if (as.table){
      networkStats <- cbind(degree(CIG), betweenness(CIG), closeness(CIG),
                            evcent(CIG)$vector, nNeg, nPos, MI, diag(S), pvars)
      colnames(networkStats) <-
        c("degree", "betweenness", "closeness", "eigenCentrality", "nNeg",
          "nPos", "mutualInfo", "variance", "partialVar")
      return(networkStats)
    }
    if (!as.table){
      return(list(degree = degree(CIG), betweenness = betweenness(CIG),
                  closeness = closeness(CIG, mode = "all"),
                  eigenCentrality = evcent(CIG, scale = FALSE)$vector,
                  nNeg = nNeg, nPos = nPos, chordal = is.chordal(CIG)$chordal,
                  mutualInfo = MI, variance = diag(S), partialVar = pvars))
    }
  }
}









#' Gaussian graphical model node pair path statistics
#' 
#' Function that calculates, for a specified node pair representing endpoints,
#' path statistics from a sparse precision matrix. The sparse precision matrix
#' is taken to represent the conditional independence graph of a Gaussian
#' graphical model. The contribution to the observed covariance between the
#' specified endpoints is calculated for each (heuristically) determined path
#' between the endpoints.
#' 
#' The conditional independence graph (as implied by the sparse precision
#' matrix) is undirected. In undirected graphs origin and destination are
#' interchangeable and are both referred to as 'endpoints' of a path. The
#' function searches for shortest paths between the specified endpoints
#' \code{node1} and \code{node2}. It searches for shortest paths that visit
#' nodes only once. The shortest paths between the provided endpoints are
#' determined heuristically by the following procedure. The search is initiated
#' by application of the \code{get.all.shortest.paths}-function from the
#' \code{\link{igraph}}-package, which yields all shortest paths between the
#' nodes. Next, the neighborhoods of the endpoints are defined (excluding the
#' endpoints themselves). Then, the shortest paths are found between: (a)
#' \code{node1} and node \emph{Vs} in its neighborhood; (b) node \emph{Vs} in
#' the \code{node1}-neighborhood and node \emph{Ve} in the
#' \code{node2}-neighborhood; and (c) node \emph{Ve} in the
#' \code{node2}-neighborhood and \code{node2}. These paths are glued and new
#' shortest path candidates are obtained (preserving only novel paths). In
#' additional iterations (specified by \code{neiExpansions}) the \code{node1}-
#' and \code{node2}-neighborhood are expanded by including their neighbors
#' (still excluding the endpoints) and shortest paths are again searched as
#' described above.
#' 
#' The contribution of a particular path to the observed covariance between the
#' specified node pair is calculated in accordance with Theorem 1 of Jones and
#' West (2005). As in Jones and West (2005), paths whose weights have an
#' opposite sign to the marginal covariance (between endnodes of the path) are
#' referred to as 'moderating paths' while paths whose weights have the same
#' sign as the marginal covariance are referred to as 'mediating' paths. Such
#' paths are visualized when \code{graph = TRUE}.
#' 
#' All arguments following the \code{graph} argument are only (potentially)
#' used when \code{graph = TRUE}. When \code{graph = TRUE} the conditional
#' independence graph is returned with the paths highlighted that have the
#' highest contribution to the marginal covariance between the specified
#' endpoints. The number of paths highlighted is indicated by \code{nrPaths}.
#' The edges of mediating paths are represented in green while the edges of
#' moderating paths are represented in red. When \code{all.edges = TRUE} the
#' edges other than those implied by the \code{nrPaths}-paths between
#' \code{node1} and node2 are also visualized (in lightgrey). When
#' \code{all.edges = FALSE} only the mediating and moderating paths implied by
#' \code{nrPaths} are visualized.
#' 
#' The default layout gives a circular placement of the vertices. Most layout
#' functions supported by \code{\link{igraph}} are supported (the function is
#' partly a wrapper around certain \code{\link{igraph}} functions). The igraph
#' layouts can be invoked by a \code{character} that mimicks a call to a
#' \code{\link{igraph}} layout functions in the \code{lay} argument. When using
#' \code{lay = NULL} one can specify the placement of vertices with the
#' \code{coords} argument. The row dimension of this matrix should equal the
#' number of (pruned) vertices. The column dimension then should equal 2 (for
#' 2D layouts) or 3 (for 3D layouts). The \code{coords} argument can also be
#' viewed as a convenience argument as it enables one, e.g., to layout a graph
#' according to the coordinates of a previous call to \code{Ugraph}. If both
#' the the lay and the coords arguments are not \code{NULL}, the lay argument
#' takes precedence
#' 
#' The arguments \code{Lcex} and \code{PTcex} are only used when \code{legend =
#' TRUE}. If \code{prune = TRUE} the vertices of degree 0 (vertices not
#' implicated by any edge) are removed. For the colors supported by the
#' arguments \code{nodecol}, \code{Vcolor}, and \code{VBcolor}, see
#' \url{https://stat.columbia.edu/~tzheng/files/Rcolor.pdf}.
#' 
#' @param P0 Sparse (possibly standardized) precision matrix.
#' @param node1 A \code{numeric} specifying an endpoint. The numeric should
#' correspond to a row/column of the precision matrix and as such represents
#' the corresponding variable.
#' @param node2 A \code{numeric} specifying a second endpoint. The numeric
#' should correspond to a row/column of the precision matrix and as such
#' represents the corresponding variable.
#' @param neiExpansions A \code{numeric} determining how many times the
#' neighborhood around the respective endpoints should be expanded in the
#' search for shortest paths between the node pair.
#' @param verbose A \code{logical} indicating if a summary of the results
#' should be printed on screen.
#' @param graph A \code{logical} indicating if the strongest paths should be
#' visualized with a graph.
#' @param nrPaths A \code{numeric} indicating the number of paths (with the
#' highest contribution to the marginal covariance between the indicated node
#' pair) to be visualized/highlighted.
#' @param lay A \code{character} mimicking a call to \code{\link{igraph}}
#' layout functions. Determines the placement of vertices.
#' @param coords A \code{matrix} containing coordinates. Alternative to the
#' lay-argument for determining the placement of vertices.
#' @param nodecol A \code{character} determining the color of \code{node1} and
#' \code{node2}.
#' @param Vsize A \code{numeric} determining the vertex size.
#' @param Vcex A \code{numeric} determining the size of the vertex labels.
#' @param VBcolor A \code{character} determining the color of the vertex
#' borders.
#' @param VLcolor A \code{character} determining the color of the vertex
#' labels.
#' @param all.edges A \code{logical} indicating if edges other than those
#' implied by the \code{nrPaths}-paths between \code{node1} and node2 should
#' also be visualized.
#' @param prune A \code{logical} determining if vertices of degree 0 should be
#' removed.
#' @param legend A \code{logical} indicating if the graph should come with a
#' legend.
#' @param scale A \code{numeric} representing a scale factor for visualizing
#' strenght of edges. It is a relative scaling factor, in the sense that the
#' edges implied by the \code{nrPaths}-paths between \code{node1} and node2
#' have edge thickness that is twice this scaling factor (so it is a scaling
#' factor vis-a-vis the unimplied edges).
#' @param Lcex A \code{numeric} determining the size of the legend box.
#' @param PTcex A \code{numeric} determining the size of the exemplary lines in
#' the legend box.
#' @param main A \code{character} giving the main figure title.
#' @return An object of class list: \item{pathStats}{A \code{matrix} specifying
#' the paths, their respective lengths, and their respective contributions to
#' the marginal covariance between the endpoints.} \item{paths}{A \code{list}
#' representing the respective paths as numeric vectors.} \item{Identifier}{A
#' \code{data.frame} in which each numeric from \code{paths} is connected to an
#' identifier such as a variable name.}
#' @note Eppstein (1998) describes a more sophisticated algorithm for finding
#' the top \emph{k} shortest paths in a graph.
#' @author Wessel N. van Wieringen, Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{ridgeP}}, \code{\link{optPenalty.LOOCVauto}},
#' \code{\link{sparsify}}
#' @references Eppstein, D. (1998). Finding the k Shortest Paths. SIAM Journal
#' on computing 28: 652-673.
#' 
#' Jones, B., and West, M. (2005). Covariance Decomposition in Undirected
#' Gaussian Graphical Models. Biometrika 92: 779-786.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p <- 25
#' n <- 10
#' set.seed(333)
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X) <- letters[1:p]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCVauto(X, lambdaMin = .5, lambdaMax = 30)
#' 
#' ## Determine support regularized standardized precision under optimal penalty
#' PC0 <- sparsify(OPT$optPrec, threshold = "localFDR")$sparseParCor
#' 
#' ## Obtain information on mediating and moderating paths between nodes 14 and 23
#' pathStats <- GGMpathStats(PC0, 14, 23, verbose = TRUE, prune = FALSE)
#' pathStats
#' 
#' @export GGMpathStats
GGMpathStats <- function(P0, node1, node2, neiExpansions = 2, verbose = TRUE,
                         graph = TRUE, nrPaths = 2, lay = "layout_in_circle",
                         coords = NULL, nodecol = "skyblue", Vsize = 15,
                         Vcex = .6, VBcolor = "darkblue", VLcolor = "black",
                         all.edges = TRUE, prune = TRUE, legend = TRUE,
                         scale = 1, Lcex = .8, PTcex = 2, main = ""){
  ##############################################################################
  # - Function that expresses the covariance between a pair of variables as a
  #   sum of path weights
  # - The sum of path weights is based on the shortest paths connecting the pair
  #   in an undirected graph
  # - P0            > sparse precision/partial correlation matrix
  # - node1         > start node of the path
  # - node2         > end node of the path
  # - neiExpansions > a numeric determining how many times the neighborhood
  #                   around the start and end node should be expanded in the
  #                   search for shortest paths between the node pair.
  #                   Default = 2
  # - verbose       > logical indicating if output should also be printed on
  #                   screen. Default = TRUE
  # - graph         > Optional argument for visualization strongest paths,
  #                   default = TRUE
  # - nrPaths	      > indicates the number of paths with the highest
  #                   contribution to the marginal covariance
  #                   between the indicated node pair (node1 and node2) to be
  #                   visualized/highlighted;
  #                   only when graph = TRUE
  # - lay           > determines layout of the graph. Most layouts in
  #                   'layout{igraph}' are accepted. Default =
  #                   layout_in_circle.
  # - coords  >       matrix of coordinates to determine layout of the graph.
  #                   The row dimension should equal the number of (pruned)
  #                   vertices. The column dimension should equal 2
  #                   (for 2D layouts) or 3 (for 3D layouts). Enables one,
  #                   e.g., to layout the graph according to the coordinates
  #                   of a previous call to Ugraph. If both the the lay and the
  #                   coords arguments are not NULL, the lay argument takes
  #                   precedence
  # - nodecol       > gives color of node1 and node2; only when graph = TRUE
  # - Vsize   	    > gives vertex size, default = 15; only when graph = TRUE
  # - Vcex    	    > gives size vertex labels, default = .6; only when
  #                   graph = TRUE
  # - VBcolor     	> gives color of the vertex border, default = "darkblue";
  #                   only when graph = TRUE
  # - VLcolor      	> gives color of the vertex labels, default = "black";
  #                   only when graph = TRUE
  # - all.edges     > logical indicating if edges other than those implied by
  #                   the 'nrPaths' paths between
  #                   node1 and node2 should also be visualized. Default = TRUE;
  #                   only when graph = TRUE
  # - prune         > logical indicating if vertices of degree 0 should be
  #                   removed. Default = TRUE; only when graph = TRUE
  # - legend        > optional inclusion of color legend, default = TRUE; only
  #                   when graph = TRUE
  # - scale         > scale factor for visualizing strenght of edges,
  #                   default = 1. It is a relative scaling
  #                   factor, in the sense that the edges implied by the
  #                   'nrPaths' paths between node1 and node2 have edge
  #                   thickness that is twice this scaling factor (so it is
  #                   a scaling factor vis-a-vis the unimplied edges); only
  #                   when all.edges = TRUE
  # - Lcex          > scaling legend box, default = .8; only when legend = TRUE
  # - PTcex         > scaling node in legend box, default = 2; only when
  #                   legend = TRUE
  # - main          > character specifying heading figure, default = ""
  #
  # - NOTES:
  # - As in Jones & West (2005), paths whose weights have an opposite sign to
  #   the marginal covariance (between endnodes of the path) are referred to
  #   as 'moderating paths' while paths whose weights have the same sign as the
  #   marginal covariance are referred to as 'mediating' paths
  ##############################################################################

  # Dependencies
  # require("base")
  # require("igraph")
  # require("reshape")

  if (!is.matrix(P0)){
    stop("Input (P0) should be a matrix")
  }
  else if (!isSymmetric(P0)){
    stop("Input (P0) should be a symmetric matrix")
  }
  else if (!evaluateS(P0, verbose = FALSE)$posEigen){
    stop("Input (P0) is expected to be positive definite")
  }
  else if (class(node1) != "numeric"){
    stop("Input (node1) is of wrong class")
  }
  else if (length(node1) != 1){
    stop("Length input (node1) must be 1")
  }
  else if (!.is.int(node1)){
    stop("Input (node1) should be a numeric integer")
  }
  else if (node1 < 1){
    stop("Input (node1) cannot be zero or negative")
  }
  else if (node1 > ncol(P0)){
    stop("Input (node1) cannot exceed the number of variables in P0")
  }
  else if (class(node2) != "numeric"){
    stop("Input (node2) is of wrong class")
  }
  else if (length(node2) != 1){
    stop("Length input (node2) must be 1")
  }
  else if (!.is.int(node2)){
    stop("Input (node2) should be a numeric integer")
  }
  else if (node2 < 1){
    stop("Input (node2) cannot be zero or negative")
  }
  else if (node2 > ncol(P0)){
    stop("Input (node2) cannot exceed the number of variables in P0")
  }
  else if (node1 == node2){
    stop("Inputs (node1 and node2) cannot be equal")
  }
  else if (class(neiExpansions) != "numeric"){
    stop("Input (neiExpansions) is of wrong class")
  }
  else if (length(neiExpansions) != 1){
    stop("Length input (neiExpansions) must be 1")
  }
  else if (!.is.int(neiExpansions)){
    stop("Input (neiExpansions) should be a numeric integer")
  }
  else if (neiExpansions < 1){
    stop("Input (neiExpansions) cannot be zero or negative")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # Some warnings
    if (all(P0 != 0)){
      warning("Given input (P0) implies a saturated conditional independence ",
              "graph")
    }
    if (all(P0[!diag(nrow(P0))] == 0)){
      warning("Given input (P0) implies an empty conditional independence ",
              "graph")
    }

    # Precision associated graph
    colnames(P0) <- colnames(P0, do.NULL = FALSE)
    Names  <- colnames(P0)
    adjMat <- adjacentMat(P0)
    colnames(adjMat) <- 1:nrow(P0)
    rownames(adjMat) <- 1:nrow(P0)
    G <- graph.adjacency(adjMat, mode = "undirected")

    # Is there a connection between the specified nodes?
    pathExits <- is.finite(shortest.paths(G, node1, node2))

    # If nodes are connected, evaluate path contributions
    if (!pathExits){
      stop("provided node pair is not connected.")
    } else {
      # Determinant of precision
      detP0 <- det(P0)

      # Objects to be returned
      paths <- list()
      pathStats <- numeric()

      # Shortest paths between the nodes
      slh <- get.all.shortest.paths(G, node1, node2)$res
      for (u in 1:length(slh)){
        fullPath <- slh[[u]]
        pName <- .path2string(fullPath)
        paths[[length(paths)+1]] <- fullPath
        pathStats <- rbind(pathStats, c(length(slh[[u]])-1,
                                        .pathContribution(P0, fullPath, detP0)))
        rownames(pathStats)[nrow(pathStats)] <- pName
      }
      nei1 <- node1
      nei2 <- node2

      for (u in 1:neiExpansions){
        # Consider longer paths between the nodes
        nei1temp <- nei1
        nei2temp <- nei2
        for (v1 in 1:length(nei1)){
          nei1temp <- c(nei1temp, neighbors(G, nei1[v1]))
        }
        for (v2 in 1:length(nei2)){
          nei2temp <- c(nei2temp, neighbors(G, nei2[v2]))
        }
        nei1 <- setdiff(unique(nei1temp), node1)
        nei2 <- setdiff(unique(nei2temp), node2)
        slh  <- .pathAndStats(G, node1, node2, nei1, nei2, P0, detP0,
                              rownames(pathStats))
        pathStats <- rbind(pathStats, slh$pathStats)
        paths <- c(paths, slh$paths)
      }

      # Wrap up
      pNames <- rownames(pathStats)
      paths  <- paths[order(abs(pathStats[,2]), decreasing=TRUE)]
      pathStats <- matrix(pathStats[order(abs(pathStats[,2]),
                                          decreasing=TRUE),], ncol = 2)
      names(paths) = rownames(pathStats) <- pNames
      colnames(pathStats) <- c("length", "contribution")
      if (verbose | graph){covNo1No2 <- solve(P0)[node1, node2]}

      # Summary
      if (verbose){
        covNo1No2expl <- sum(pathStats[,2])

        # Reformat results
        statsTable <- data.frame(cbind(rownames(pathStats), pathStats[,1],
                                       round(pathStats[,2], 5)))
        rownames(statsTable) <- NULL
        colnames(statsTable) <- c("path", "length", "contribution")

        # Print results on screen
        cat("Covariance between node pair :", round(covNo1No2, 5), "\n")
        cat("----------------------------------------\n")
        print(statsTable, quote=FALSE)
        cat("----------------------------------------\n")
        cat("Sum path contributions       :", round(covNo1No2expl, 5), "\n")

      }

      # Visualize
      if (graph){
        if (class(nrPaths) != "numeric"){
          stop("Input (nrPaths) is of wrong class")
        } else if (length(nrPaths) != 1){
          stop("Length input (nrPaths) must be one")
        } else if (!.is.int(nrPaths)){
          stop("Input (nrPaths) should be a numeric integer")
        } else if (nrPaths <= 0){
          stop("Input (nrPaths) must be a strictly positive integer")
        } else if (nrPaths > length(paths)){
          stop("Input (nrPaths) cannot exceed the total number of paths ",
               "discerned")
        } else if (!((length(intersect(lay,c("layout_as_star", "layout_as_tree",
                                             "layout_in_circle", "layout_nicely",
                                             "layout_with_dh", "layout_with_gem",
                                             "layout_with_graphopt", "layout_on_grid",
                                             "layout_with_mds", "layout_components",
                                             "layout_on_sphere", "layout_randomly",
                                             "layout_with_fr", "layout_with_kk",
                                             "layout_with_lgl"))) > 0) | is.null(lay))){
          stop("lay should be 'NULL' or one of
               {'layout_as_star', 'layout_as_tree',
               'layout_in_circle', 'layout_nicely',
               'layout_with_dh', 'layout_with_gem',
               'layout_with_graphopt', 'layout_on_grid',
               'layout_with_mds', 'layout_components',
               'layout_on_sphere', 'layout_randomly',
               'layout_with_fr', 'layout_with_kk',
               'layout_with_lgl'}")
        } else if (!is.null(coords) & !inherits(coords, "matrix")){
          stop("Input (coords) is of wrong class")
      } else if (is.null(lay) & is.null(coords)){
        stop("Input (lay) and input (coords) cannot be both NULL")
      } else if (class(nodecol) != "character"){
        stop("Input (nodecol) is of wrong class")
      } else if (length(nodecol) != 1){
        stop("Length input (nodecol) must be one")
      } else if (class(Vsize) != "numeric"){
        stop("Input (Vsize) is of wrong class")
      } else if (length(Vsize) != 1){
        stop("Length input (Vsize) must be one")
      } else if (Vsize <= 0){
        stop("Input (Vsize) must be strictly positive")
      } else if (class(Vcex) != "numeric"){
        stop("Input (Vcex) is of wrong class")
      } else if (length(Vcex) != 1){
        stop("Length input (Vcex) must be one")
      } else if (Vcex <= 0){
        stop("Input (Vcex) must be strictly positive")
      } else if (class(VBcolor) != "character"){
        stop("Input (VBcolor) is of wrong class")
      } else if (length(VBcolor) != 1){
        stop("Length input (VBcolor) must be one")
      } else if (class(VLcolor) != "character"){
        stop("Input (VLcolor) is of wrong class")
      } else if (length(VLcolor) != 1){
        stop("Length input (VLcolor) must be one")
      } else if (class(all.edges) != "logical"){
        stop("Input (all.edges) is of wrong class")
      } else if (class(prune) != "logical"){
        stop("Input (prune) is of wrong class")
      } else if (class(legend) != "logical"){
        stop("Input (legend) is of wrong class")
      } else if (class(main) != "character"){
        stop("Input (main) is of wrong class")
      } else {
        # Preliminaries
        AM <- adjacentMat(P0)
        GA <- graph.adjacency(AM, mode = "undirected")
        if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}
        colnames(P0) = rownames(P0) <- seq(1, ncol(P0), by = 1)
        Mmelt <- melt(P0)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]

        # Layout specification
        if(is.null(lay)){
          if(dim(coords)[1] != length(V(GA))){
            stop("Row dimension of input (coords) does not match the
                 number of vertices to be plotted")
          } else if (dim(coords)[2] > 3){
            stop("Column dimension of input (coords) exceeds the number
                 of dimensions that can be visualized")
          } else {lays = coords}
          }
        else{
          if(lay == "layout_as_star"){
            lays = igraph::layout_as_star(GA)}
          if(lay == "layout_as_tree")
          {lays = igraph::layout_as_tree(GA)}
          if(lay == "layout_in_circle"){
            lays = igraph::layout_in_circle(GA)}
          if(lay == "layout_nicely"){
            lays = igraph::layout_nicely(GA)}
          if(lay == "layout_with_dh"){
            lays = igraph::layout_with_dh(GA)}
          if(lay == "layout_with_gem"){
            lays = igraph::layout_with_gem(GA)}
          if(lay == "layout_with_graphopt"){
            lays = igraph::layout_with_graphopt(GA)}
          if(lay == "layout_on_grid"){
            lays = igraph::layout_on_grid(GA)}
          if(lay == "layout_with_mds"){
            lays = igraph::layout_with_mds(GA)}
          if(lay == "layout_components"){
            lays = igraph::layout_components(GA)}
          if(lay == "layout_on_sphere"){
            lays = igraph::layout_on_sphere(GA)}
          if(lay == "layout_randomly"){
            lays = igraph::layout_randomly(GA)}
          if(lay == "layout_with_fr"){
            lays = igraph::layout_with_fr(GA)}
          if(lay == "layout_with_kk"){
            lays = igraph::layout_with_kk(GA)}
          if(lay == "layout_with_lgl"){
            lays = igraph::layout_with_lgl(GA)}
        }

        # Determine if path is mediating or moderating and color accordingly
        for (k in 1:nrPaths){
          Path <- unlist(paths[k])
          if (sign(covNo1No2) == sign(pathStats[k,2])){COL = "green"}
          if (sign(covNo1No2) != sign(pathStats[k,2])){COL = "red"}
          for (i in 1:(length(Path) - 1)){
            if (Path[i] > Path[i + 1]){
              tempX1 <- Path[i]
              tempX2 <- Path[i + 1]
            } else {
              tempX1 <- Path[i + 1]
              tempX2 <- Path[i]
            }
            row <- which(Mmelt$X1 == tempX1 & Mmelt$X2 == tempX2)
            E(GA)[row]$color <- COL
          }
        }

        # Coloring nodes
        V(GA)$color <- "white"
        V(GA)$color[node1] <- nodecol
        V(GA)$color[node2] <- nodecol

        # Produce graph
        if (all.edges){
          if (class(scale) != "numeric"){
            stop("Input (scale) is of wrong class")
          } else if (length(scale) != 1){
            stop("Length input (scale) must be one")
          } else if (scale <= 0){
            stop("Input (scale) must be strictly positive")
          } else {
            E(GA)[is.na(E(GA)$color)]$color <- "grey"
            E(GA)$weight <- 1
            E(GA)[E(GA)$color == "green"]$weight <- 2
            E(GA)[E(GA)$color == "red"]$weight   <- 2
            plot(GA, layout = lays, vertex.size = Vsize,
                 vertex.label.family = "sans", vertex.label.cex = Vcex,
                 edge.width = scale*abs(E(GA)$weight),
                 vertex.color = V(GA)$color, vertex.frame.color = VBcolor,
                 vertex.label.color = VLcolor, main = main)
          }
        } else {
          plot(GA, layout = lays, vertex.size = Vsize,
               vertex.label.family = "sans", vertex.label.cex = Vcex,
               vertex.color = V(GA)$color, vertex.frame.color = VBcolor,
               vertex.label.color = VLcolor, main = main)
        }

        # Legend
        if (legend){
          if (class(Lcex) != "numeric"){
            stop("Input (Lcex) is of wrong class")
          } else if (length(Lcex) != 1){
            stop("Length input (Lcex) must be one")
          } else if (Lcex <= 0){
            stop("Input (Lcex) must be strictly positive")
          } else if (class(PTcex) != "numeric"){
            stop("Input (PTcex) is of wrong class")
          } else if (length(PTcex) != 1){
            stop("Length input (PTcex) must be one")
          } else if (PTcex <= 0){
            stop("Input (PTcex) must be strictly positive")
          } else {
            legend("bottomright", c("mediating path", "moderating path"),
                   lty=c(1,1), col = c("green", "red"), cex = Lcex,
                   pt.cex = PTcex)
          }
        }
      }
    }

    # Return
    Numeric    <- rownames(adjMat)
    VarName    <- Names
    identifier <- data.frame(Numeric, VarName)
    return(list(pathStats = pathStats, paths = paths,
                  Identifier = identifier))
    }
  }
}




##------------------------------------------------------------------------------
##
## Wrapper function
##
##------------------------------------------------------------------------------







#' Wrapper function
#' 
#' Function that forms a wrapper around certain \code{rags2ridges}
#' functionalities. More specifically, it (automatically) invokes
#' functionalities to get from high-dimensional data to a penalized precision
#' estimate, to the corresponding conditional independence graph and topology
#' summaries.
#' 
#' The wrapper always uses the alternative ridge precision estimator (see
#' \code{\link{ridgeP}}) with \code{target} as the target matrix. The optimal
#' value for the penalty parameter is determined by employing Brent's method to
#' the calculation of a cross-validated negative log-likelihood score (see
#' \code{\link{optPenalty.LOOCVauto}}). The support of the regularized
#' precision matrix is determined by way of local FDR thresholding (see
#' \code{\link{sparsify}}). The corresponding conditional independence graph is
#' visualized using \code{\link{Ugraph}} with \code{type = "fancy"}. This
#' visualization as well as the calculation of network statistics (see
#' \code{\link{GGMnetworkStats}}) is based on the standardization of the
#' regularized and sparsified precision matrix to a partial correlation matrix.
#' 
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param dir A \code{character} specifying the directory in which the (visual)
#' output is to be stored.
#' @param fileTypeFig A \code{character} determining the file type of visual
#' output. Must be one of: "pdf", "eps".
#' @param FDRcut A \code{numeric} indicating the cut-off for partial
#' correlation element selection based on local FDR thresholding.
#' @param nOutput A \code{logical} indicating if numeric output should be
#' returned.
#' @param verbose A \code{logical} indicating if progress updates should be
#' printed on screen.
#' @return The function stores in the specified directory \code{dir} a
#' condition number plot (either .pdf or .eps file), a visualization of the
#' network (either .pdf or .eps file), and a file containing network statistics
#' (.txt file). When \code{nOutput = TRUE} the function also returns an object
#' of class \code{list}: \item{optLambda}{A \code{numeric} giving the optimal
#' value of the penalty parameter.} \item{optPrec}{A \code{matrix} representing
#' the regularized precision matrix under the optimal value of the penalty
#' parameter.} \item{sparseParCor}{A \code{matrix} representing the sparsified
#' partial correlation matrix.} \item{networkStats}{A \code{matrix} giving the
#' calculated network statistics.}
#' @note We consider this to be a preliminary version of an envisioned wrapper
#' than will take better form with subsequent versions of \code{rags2ridges}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N. van Wieringen
#' @seealso \code{\link{ridgeP}}, \code{\link{conditionNumberPlot}},
#' \code{\link{optPenalty.LOOCVauto}}, \code{\link{sparsify}},
#' \code{\link{Ugraph}}, \code{\link{GGMnetworkStats}}
#' @examples
#' 
#' \dontrun{
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Employ the wrapper function
#' theWorks <- fullMontyS(X, lambdaMin = .5, lambdaMax = 30)}
#' 
#' @export fullMontyS
fullMontyS <- function(Y, lambdaMin, lambdaMax,
                       target = default.target(covML(Y)), dir = getwd(),
                       fileTypeFig = "pdf", FDRcut = .9, nOutput = TRUE,
                       verbose = TRUE){
  ##############################################################################
  # - Function that forms a wrapper around the rags2ridges functionalities
  # - Invokes functionalities to get from data to graph and topology summaries
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter
  # - lambdaMax   > maximum value penalty parameter
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - dir         > specifies the directory in which the (visual) output is
  #                 stored
  # - fileTypeFig > signifies filetype of visual output; Should be one of
  #                 {"pdf", "eps"}
  # - FDRcut      > cut-off for partial correlation element selection based on
  #                 local FDR thresholding. Default = .9.
  # - nOutput     > logical indicating if numeric output should be given
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  #
  # - NOTES:
  # - Always uses the alternative ridge estimator by van Wieringen and Peeters
  #   (2015)
  # - Always uses LOOCV by Brent for optimal penalty parameter determination
  # - Always uses support determination by local FDR thresholding
  # - Visualizes the network by 'Ugraph' with type = "fancy". Network on basis
  #   partial correlations
  # - Network statistics calculated on sparsified partial correlation network
  # - There are no elaborate input checks as these are all covered by the
  #   indvidual functions invoked
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("graphics")
  # require("Hmisc")
  # require("fdrtool")
  # require("igraph")
  # require("reshape")
  # require("utils")

  if (class(dir) != "character"){
    stop("Specify directory for output (dir) as 'character'")
  }
  else if (!(fileTypeFig %in% c("pdf", "eps"))){
    stop("Input (fileTypeFig) should be one of {'pdf', 'eps'}")
  }
  else if (class(nOutput) != "logical"){
    stop("Input (nOutput) is of wrong class")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # In case one does not know:
    cat("Output files are stored in the following directory:", dir, "\n")
    cat("\n")

    # Determine optimal penalty and precision
    if (verbose){cat("Progress:", "\n")}
    if (verbose){cat("Determining optimal penalty value...", "\n")}
    optimal <- optPenalty.LOOCVauto(Y, lambdaMin = lambdaMin,
                                    lambdaMax = lambdaMax,
                                    lambdaInit = (lambdaMin+lambdaMax)/2,
                                    target = target)

    # Condition number plot
    if (verbose){cat("Generating condition number plot...", "\n")}
    if (fileTypeFig == "pdf"){pdf(paste(dir, "Condition_Number_Plot.pdf"))}
    if (fileTypeFig == "eps"){
      setEPS()
      postscript(paste(dir, "Condition_Number_Plot.eps"))
    }
    conditionNumberPlot(covML(Y), lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                        step = 100000, target = target, vertical = TRUE,
                        value = optimal$optLambda, main = FALSE, verbose =FALSE)
    dev.off()

    # Sparsify the precision matrix
    if (verbose){cat("Determining support...", "\n")}
    PC0 <- sparsify(optimal$optPrec, threshold = "localFDR", FDRcut = FDRcut,
                    output = "heavy", verbose = FALSE)$sparseParCor

    # Visualize the network
    if (verbose){cat("Visualizing network...", "\n")}
    if (fileTypeFig == "pdf"){pdf(paste(dir, "Network.pdf"))}
    if (fileTypeFig == "eps"){setEPS(); postscript(paste(dir, "Network.eps"))}
    Ugraph(PC0, type = "fancy", Vsize = ncol(PC0)/10, Vcex = ncol(PC0)/130)
    dev.off()

    # Calculate network statistics
    if (verbose){cat("Calculating network statistics...", "\n")}
    Stats <- GGMnetworkStats(PC0, as.table = TRUE)
    capture.output(print(Stats), file = paste(dir, "Network_Statistics.txt"))

    # Done
    if (verbose){cat("DONE!", "\n")}

    # Return
    if (nOutput){
      return(list(optLambda = optimal$optLambda, optPrec = optimal$optPrec,
                  sparseParCor = PC0, networkStats = Stats))
    }
  }
}






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module B: rags2ridges Fused
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# See R/rags2ridgesFused.R






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module C: rags2ridges Miscellaneous
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# See R/rags2ridgesMisc.R








#' Moments of the sample covariance matrix.
#' 
#' Calculates the moments of the sample covariance matrix. It assumes that the
#' summands (the outer products of the samples' random data vector) that
#' constitute the sample covariance matrix follow a Wishart-distribution with
#' scale parameter \eqn{\mathbf{\Sigma}} and shape parameter \eqn{\nu}. The
#' latter is equal to the number of summands in the sample covariance estimate.
#' 
#' 
#' @param Sigma Positive-definite \code{matrix}, the scale parameter
#' \eqn{\mathbf{\Sigma}} of the Wishart distribution.
#' @param shape A \code{numeric}, the shape parameter \eqn{\nu} of the Wishart
#' distribution. Should exceed the number of variates (number of rows or
#' columns of \code{Sigma}).
#' @param moment An \code{integer}. Should be in the set \eqn{\{-4, -3, -2, -1,
#' 0, 1, 2, 3, 4\}} (only those are explicitly specified in Lesac, Massam,
#' 2004).
#' @return The \eqn{r}-th moment of a sample covariance matrix:
#' \eqn{E(\mathbf{S}^r)}.
#' @author Wessel N. van Wieringen.
#' @references Lesac, G., Massam, H. (2004), "All invariant moments of the
#' Wishart distribution", \emph{Scandinavian Journal of Statistics}, 31(2),
#' 295-318.
#' @examples
#' 
#' # create scale parameter
#' Sigma <- matrix(c(1, 0.5, 0, 0.5, 1, 0, 0, 0, 1), byrow=TRUE, ncol=3)
#' 
#' # evaluate expectation of the square of a sample covariance matrix 
#' # that is assumed to Wishart-distributed random variable with the 
#' # above scale parameter Sigma and shape parameter equal to 40.
#' momentS(Sigma, 40, 2)
#' 
#' @export momentS
momentS <- function(Sigma,
                    shape,
                    moment=1){

	########################################################################
	#
	# DESCRIPTION:
	# Returns the moments of a Wishart-distributed random variable.
	# Only those explicitly given in Lesac, Massam (2004) are implemented.
	#
	# ARGUMENTS:
	# -> Sigma	: Positive-definite 'matrix', the scale parameter of
	#                 the Wishart distribution.
	# -> shape	: A 'numeric', the shape parameter of the Wishart
	#                 distribution. Should exceed the number of variates.
	# -> moment	: An 'integer'. Should be in the set
	#                 {-4, -3, -2, -1, 0, 1, 2, 3, 4} (only those are
	#                 explicitly specified in Lesac, Massam, 2004).
	#
	# DEPENDENCIES:
	# Currently, none.
	#
	# NOTES:
   	# ....
	#
	########################################################################

	# input checks
	if (shape < nrow(Sigma) + 1){
		stop("shape parameter should exceed the number of variates")
	}
	if (all(moment != c(-c(4:1), 0:4))){
		stop("moment not implemented")
	}

	# number of variates
	p <- nrow(Sigma)

	# moment, case-wise
	if (moment==-4){
		Sinv     <- solve(Sigma);
		constant <- (shape-p+2) * (shape-p+1) * (shape-p-2) * (shape-p) *
		            (shape-p-7) * (shape-p-5) * (shape-p-3) * (shape-p-1);
		c1       <- 5 * (shape-p) - 11;
		c2       <- 5 * (shape-p)^2 - 16*(shape-p) + 11;
		c3       <- (shape-p)^3 - 4*(shape-p)^2 + 7*(shape-p) - 4;
		c4       <- 2*(shape-p)^3 - 9*(shape-p)^2  + 24*(shape-p) + 19;
		c5       <- (shape-p)^4 - 5*(shape-p)^3 +
		            11*(shape-p)^2 - 11*(shape-p) + 4;
		ESr      <- (c1 * Sinv * (sum(diag(Sinv)))^3  +
		             c2 * (Sinv * sum(diag(Sinv)) * sum(Sinv * Sinv) +
		             Sinv %*% Sinv * (sum(diag(Sinv)))^2) +
		             c3 * (Sinv * sum(diag(Sinv %*% Sinv %*% Sinv)) +
		             3 * Sinv %*% Sinv %*% Sinv * sum(diag(Sinv))) +
		             c4 * Sinv %*% Sinv * sum(Sinv * Sinv) +
		             c5 * Sinv %*% Sinv %*% Sinv %*% Sinv) / constant
	}
	if (moment==-3){
		Sinv     <- solve(Sigma);
		constant <- ((shape - p) * (shape - p - 1) * (shape - p - 3) *
		             (shape - p + 1) * (shape - p - 5))
		ESr      <- ((2 * (.trace(Sinv))^2 * Sinv +
		             (shape - p - 1) * (.trace(Sinv %*% Sinv) * Sinv +
		             2 * .trace(Sinv) * Sinv %*% Sinv) +
		             (shape - p - 1)^2 * Sinv %*% Sinv %*% Sinv) / constant)
	}
	if (moment==-2){
		Sinv <- solve(Sigma);
		ESr  <- (((shape - p - 1) * Sinv %*% Sinv
		          + .trace(Sinv) * Sinv) /
			  ((shape - p) * (shape - p - 1) * (shape - p - 3)))
	}
	if (moment==-1){
		ESr <- solve(Sigma) / (shape - p - 1)
	}
	if (moment==0){
		ESr <- diag(p)
	}
	if (moment==1){
		ESr <- shape * Sigma
	}
	if (moment==2){
		ESr <- (shape * (shape + 1) * Sigma %*% Sigma +
			shape * .trace(Sigma) * Sigma)
	}
	if (moment==3){
		ESr <- (shape * (.trace(Sigma))^2 * Sigma +
			shape * (shape + 1) * (.trace(Sigma %*% Sigma) * Sigma +
			2 * .trace(Sigma) * Sigma %*% Sigma) +
			shape * (shape^2 + 3*shape + 4) * Sigma %*% Sigma %*% Sigma)
	}
	if (moment==4){
		ESr <- (shape * .trace(Sigma %*% Sigma %*% Sigma) * Sigma +
			3 * shape * (shape + 1) * (.trace(Sigma) * .trace(Sigma %*% Sigma) * Sigma +
                            (.trace(Sigma))^2 * Sigma %*% Sigma) +
                            shape * (shape^2 + 3 * shape + 4) * (3* .trace(Sigma) * Sigma %*% Sigma %*% Sigma +
			    .trace(Sigma %*% Sigma %*% Sigma) * Sigma) +
                            shape * (2 * shape^2 + 5 * shape + 5) * .trace(Sigma %*% Sigma) * Sigma %*% Sigma +
                            shape * (shape^3 + 6 * shape^2 + 21 * shape + 20) * Sigma %*% Sigma %*% Sigma %*% Sigma)
	}
	return(ESr / shape^moment)
}









#' Prune square matrix to those variables having nonzero entries
#' 
#' Convenience function that prunes a square matrix to those variables
#' (features) having nonzero row (column) entries (i.e., to features implied in
#' graphical connections).
#' 
#' 
#' @param M (Possibly sparsified) square \code{matrix}.
#' @return A pruned \code{matrix}.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)
#' 
#' ## Determine support regularized standardized precision under optimal penalty
#' PC0 <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor
#' 
#' ## Prune sparsified partial correlation matrix
#' PC0P <- pruneMatrix(PC0)
#' 
#' @export pruneMatrix
pruneMatrix <- function(M){
  ##############################################################################
  # - Function that prunes a matrix to those variables implied in edges
  # - M > (Possibly sparsified) precision matrix
  ##############################################################################

  # Dependencies
  # require("igraph")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be square matrix")
  }
  else {
    # Preliminaries
    AM <- adjacentMat(M)
    GA <- graph.adjacency(AM, mode = "undirected")
    GA <- delete.vertices(GA, which(degree(GA) < 1))

    # Prune
    Mprune <- M[rownames(M) %in% V(GA)$name, colnames(M) %in% V(GA)$name]

    # Return
    return(Mprune)
  }
}









#' Subset 2 square matrices to union of variables having nonzero entries
#' 
#' Convenience function that subsets 2 square matrices (over the same features)
#' to the union of features that have nonzero row (column) entries (i.e.,
#' features implied in graphical connections).
#' 
#' Say you have 2 class-specific precision matrices that are estimated over the
#' same variables/features. For various reasons (such as, e.g., the desire to
#' visualize pruned class-specific networks in the same coordinates) one may
#' want to prune these matrices to those features that are implied in graphical
#' connections in at least 1 class.
#' 
#' @param M1 (Possibly sparsified) square \code{matrix}.
#' @param M2 (Possibly sparsified) square \code{matrix} over the same features
#' as \code{M1}.
#' @return An object of class list: \item{M1subset}{A pruned \code{matrix} for
#' class 1.} \item{M2subset}{A pruned \code{matrix} for class 2.}
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{Ugraph}}
#' @examples
#' 
#' ## Invoke data
#' data(ADdata)
#' 
#' ## Subset
#' ADclass1 <- ADmetabolites[, sampleInfo$ApoEClass == "Class 1"]
#' ADclass2 <- ADmetabolites[, sampleInfo$ApoEClass == "Class 2"]
#' 
#' ## Transpose data
#' ADclass1 <- t(ADclass1)
#' ADclass2 <- t(ADclass2)
#' 
#' ## Correlations for subsets
#' rAD1 <- cor(ADclass1)
#' rAD2 <- cor(ADclass2)
#' 
#' ## Simple precision estimates
#' P1 <- ridgeP(rAD1, 2)
#' P2 <- ridgeP(rAD2, 2)
#' Plist = list(P1 = P1, P2 = P2)
#' 
#' ## Threshold matrices
#' Mats <- sparsify.fused(Plist, threshold = "top", top = 20)
#' 
#' ## Prune sparsified partial correlation matrices
#' ## To union of features implied by edge
#' MatsPrune <- Union(Mats$P1$sparseParCor, Mats$P2$sparseParCor)
#' 
#' @export Union
Union <- function(M1, M2){
  ##############################################################################
  # - Function that subsets square matrices to union of features implied
  #   in edges
  # - M1 > Sparsified (precision) matrix
  # - M2 > Sparsified (precision) matrix
  ##############################################################################

  # Dependencies
  # require("reshape")

  if (!is.matrix(M1)){
    stop("M1 should be a matrix")
  }
  else if (nrow(M1) != ncol(M1)){
    stop("M1 should be square matrix")
  }
  else if (!is.matrix(M2)){
    stop("M2 should be a matrix")
  }
  else if (nrow(M2) != ncol(M2)){
    stop("M2 should be square matrix")
  }
  else {
    # Unique features matrix 1
    AmatM1  <- adjacentMat(M1)
    Mmelt   <- melt(AmatM1)
    Mmelt   <- Mmelt[Mmelt$value != 0,]
    whichM1 <- unique(c(as.character(Mmelt$X1),as.character(Mmelt$X2)))

    # Unique features matrix 2
    AmatM2  <- adjacentMat(M2)
    Mmelt   <- melt(AmatM2)
    Mmelt   <- Mmelt[Mmelt$value != 0,]
    whichM2 <- unique(c(as.character(Mmelt$X1),as.character(Mmelt$X2)))

    # Union of features present (in terms of implied in edge) in all matrices
    Union <- union(whichM1, whichM2)

    # Subset matrices
    M1s <- M1[rownames(M1) %in% Union, colnames(M1) %in% Union]
    M2s <- M2[rownames(M2) %in% Union, colnames(M2) %in% Union]

    # Return
    return(list(M1subset = M1s, M2subset = M2s))
  }
}









#' Search and visualize community-structures
#' 
#' Function that searches for and visualizes community-structures in graphs.
#' 
#' Communities in a network are groups of vertices (modules) that are densely
#' connected within. Community search is performed by the Girvan-Newman
#' algorithm (Newman and Girvan, 2004).
#' 
#' When \code{graph = TRUE} the community structure in the graph is visualized.
#' The default layout is according to the Fruchterman-Reingold algorithm
#' (1991). Most layout functions supported by \code{\link{igraph}} are
#' supported (the function is partly a wrapper around certain
#' \code{\link{igraph}} functions). The igraph layouts can be invoked by a
#' \code{character} that mimicks a call to a \code{\link{igraph}} layout
#' functions in the \code{lay} argument. When using \code{lay = NULL} one can
#' specify the placement of vertices with the \code{coords} argument. The row
#' dimension of this matrix should equal the number of vertices. The column
#' dimension then should equal 2 (for 2D layouts) or 3 (for 3D layouts). The
#' \code{coords} argument can also be viewed as a convenience argument as it
#' enables one, e.g., to layout a graph according to the coordinates of a
#' previous call to \code{Ugraph}. If both the the lay and the coords arguments
#' are not \code{NULL}, the lay argument takes precedence. Communities are
#' indicated by color markings.
#' 
#' @param P Sparsified precision \code{matrix}
#' @param graph A \code{logical} indicating if the results should be
#' visualized.
#' @param lay A \code{character} mimicking a call to \code{\link{igraph}}
#' layout functions. Determines the placement of vertices.
#' @param coords A \code{matrix} containing coordinates. Alternative to the
#' lay-argument for determining the placement of vertices.
#' @param Vsize A \code{numeric} determining the vertex size.
#' @param Vcex A \code{numeric} determining the size of the vertex labels.
#' @param Vcolor A \code{character} (scalar or vector) determining the vertex
#' color.
#' @param VBcolor A \code{character} determining the color of the vertex
#' border.
#' @param VLcolor A \code{character} determining the color of the vertex
#' labels.
#' @param main A \code{character} giving the main figure title.
#' @return An object of class list: \item{membership}{\code{numeric} vector
#' indicating, for each vertex, community membership.}
#' \item{modularityscore}{\code{numeric} scalar indicating the modularity value
#' of the community structure.}
#' 
#' When \code{graph = TRUE} the function also returns a graph.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{Ugraph}}
#' @references Csardi, G. and Nepusz, T. (2006). The igraph software package
#' for complex network research. InterJournal, Complex Systems 1695.
#' http://igraph.sf.net
#' 
#' Fruchterman, T.M.J., and Reingold, E.M. (1991). Graph Drawing by
#' Force-Directed Placement. Software: Practice & Experience, 21: 1129-1164.
#' 
#' Newman, M. and Girvan, M. (2004). Finding and evaluating community structure
#' in networks. Physical Review E, 69: 026113.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty
#' OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)
#' 
#' ## Determine support regularized standardized precision under optimal penalty
#' PC0 <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor
#' 
#' ## Search and visualize communities
#' Commy <- Communities(PC0)
#' 
#' @export Communities
Communities <- function(P, graph = TRUE, lay = "layout_with_fr", coords = NULL,
                        Vsize = 15, Vcex = 1, Vcolor = "orangered",
                        VBcolor = "darkred", VLcolor = "black", main = ""){
  ##############################################################################
  # - Function that computes and visualizes community-structures
  # - Function is partly a wrapper around certain 'igraph' functions
  # - P       > Sparsified precision matrix
  # - graph   > Logical indicating if also a graph should be returned.
  # - lay     > determines layout of the graph. Most layouts in 'layout{igraph}'
  #             are accepted. Default = layout_with_fr.
  #             Only used when graph = TRUE.
  # - coords  > matrix of coordinates to determine layout of the graph.
  #             The row dimension should equal the number of (pruned) vertices.
  #             The column dimension should equal 2 (for 2D layouts) or
  #             3 (for 3D layouts). Enables one, e.g., to layout the graph
  #             according to the coordinates of a previous call to Ugraph.
  #             If both the the lay and the coords arguments are not NULL,
  #             the lay argument takes precedence
  #             Only used when graph = TRUE.
  # - Vsize   > gives vertex size, default = 15
  #             Only used when graph = TRUE.
  # - Vcex    > gives size vertex labels, default = 1
  #             Only used when graph = TRUE.
  # - Vcolor  > gives vertex color, default = "orangered", must be character.
  #             May also be a character vector.
  #             Only used when graph = TRUE.
  # - VBcolor > gives color of the vertex border, default = "darkred"
  #             Only used when graph = TRUE.
  # - VLcolor > gives color of the vertex labels, default = "black"
  #             Only used when graph = TRUE.
  # - main    > character specifying heading figure, default = ""
  #             Only used when graph = TRUE.
  #
  # NOTES
  # - Communities on the basis of the assumption of undirected graphs
  # - Will be expanded and bettered
  ##############################################################################

  # Dependencies
  # require("base")
  # require("igraph")
  # require("reshape")

  if (!is.matrix(P)){
    stop("Input (P) should be a matrix")
  }
  else if (nrow(P) != ncol(P)){
    stop("Input (P) should be square matrix")
  }
  else {
    # Preliminaries
    AM <- adjacentMat(P)
    GA <- graph.adjacency(AM, mode = "undirected")

    # Find communities
    Communities <- cluster_edge_betweenness(GA, weights = NULL, directed = FALSE)

    # Visualize
    if(graph){
      if (!((length(intersect(lay,c("layout_as_star", "layout_as_tree",
                                    "layout_in_circle", "layout_nicely",
                                    "layout_with_dh", "layout_with_gem",
                                    "layout_with_graphopt", "layout_on_grid",
                                    "layout_with_mds", "layout_components",
                                    "layout_on_sphere", "layout_randomly",
                                    "layout_with_fr", "layout_with_kk",
                                    "layout_with_lgl"))) > 0) | is.null(lay))){
        stop("lay should be 'NULL' or one of
             {'layout_as_star', 'layout_as_tree',
             'layout_in_circle', 'layout_nicely',
             'layout_with_dh', 'layout_with_gem',
             'layout_with_graphopt', 'layout_on_grid',
             'layout_with_mds', 'layout_components',
             'layout_on_sphere', 'layout_randomly',
             'layout_with_fr', 'layout_with_kk',
             'layout_with_lgl'}")
      }
      else if (!is.null(coords) & !inherits(coords, "matrix")){
        stop("Input (coords) is of wrong class")
      }
      else if (is.null(lay) & is.null(coords)){
        stop("Input (lay) and input (coords) cannot be both NULL")
      }
      else if (class(Vsize) != "numeric"){
        stop("Input (Vsize) is of wrong class")
      }
      else if (length(Vsize) != 1){
        stop("Length Vsize must be one")
      }
      else if (Vsize <= 0){
        stop("Vsize must be positive")
      }
      else if (class(Vcex) != "numeric"){
        stop("Input (Vcex) is of wrong class")
      }
      else if (length(Vcex) != 1){
        stop("Length Vcex must be one")
      }
      else if (Vcex <= 0){
        stop("Vcex must be positive")
      }
      else if (class(Vcolor) != "character"){
        stop("Input (Vcolor) is of wrong class")
      }
      else if (length(Vcolor) != 1 & length(Vcolor) != nrow(P)){
        stop("Length Vcolor must be either one
             or equal to row (or column) dimension of P")
      }
      else if (class(VBcolor) != "character"){
        stop("Input (VBcolor) is of wrong class")
      }
      else if (length(VBcolor) != 1){
        stop("Length VBcolor must be one")
      }
      else if (class(VLcolor) != "character"){
        stop("Input (VLcolor) is of wrong class")
      }
      else if (length(VLcolor) != 1){
        stop("Length VLcolor must be one")
      }
      else if (class(main) != "character"){
        stop("Input (main) is of wrong class")
      }
      else {
        # Layout specification
        if(is.null(lay)){
          if(dim(coords)[1] != length(V(GA))){
            stop("Row dimension of input (coords) does not match the
                 number of vertices to be plotted")
          } else if (dim(coords)[2] > 3){
            stop("Column dimension of input (coords) exceeds the number
                 of dimensions that can be visualized")
          } else {lays = coords}
          }
        else{
          if(lay == "layout_as_star"){
            lays = igraph::layout_as_star(GA)}
          if(lay == "layout_as_tree")
          {lays = igraph::layout_as_tree(GA)}
          if(lay == "layout_in_circle"){
            lays = igraph::layout_in_circle(GA)}
          if(lay == "layout_nicely"){
            lays = igraph::layout_nicely(GA)}
          if(lay == "layout_with_dh"){
            lays = igraph::layout_with_dh(GA)}
          if(lay == "layout_with_gem"){
            lays = igraph::layout_with_gem(GA)}
          if(lay == "layout_with_graphopt"){
            lays = igraph::layout_with_graphopt(GA)}
          if(lay == "layout_on_grid"){
            lays = igraph::layout_on_grid(GA)}
          if(lay == "layout_with_mds"){
            lays = igraph::layout_with_mds(GA)}
          if(lay == "layout_components"){
            lays = igraph::layout_components(GA)}
          if(lay == "layout_on_sphere"){
            lays = igraph::layout_on_sphere(GA)}
          if(lay == "layout_randomly"){
            lays = igraph::layout_randomly(GA)}
          if(lay == "layout_with_fr"){
            lays = igraph::layout_with_fr(GA)}
          if(lay == "layout_with_kk"){
            lays = igraph::layout_with_kk(GA)}
          if(lay == "layout_with_lgl"){
            lays = igraph::layout_with_lgl(GA)}
        }

        # Plot
        Names <- colnames(P)
        colnames(P) = rownames(P) <- seq(1, ncol(P), by = 1)
        Mmelt <- melt(P)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]
        E(GA)$weight <- Mmelt$value
        E(GA)$color  <- "black"
        E(GA)[E(GA)$weight < 0]$style <- "dashed"
        E(GA)[E(GA)$weight > 0]$style <- "solid"
        plot(Communities, GA,
             layout = lays,
             vertex.size = Vsize,
             vertex.label.family = "sans",
             vertex.label.cex = Vcex,
             vertex.color = Vcolor,
             vertex.frame.color = VBcolor,
             vertex.label.color = VLcolor,
             edge.color = E(GA)$color,
             edge.lty = E(GA)$style,
             main = main)
      }
    }

    # Return
    return(list(membership = membership(Communities),
                modularityscore = modularity(Communities)))
  }
}









#' Visualize the differential graph
#' 
#' Function visualizing the differential graph, i.e., the network of edges that
#' are unique for 2 class-specific graphs over the same vertices
#' 
#' Say you have 2 class-specific precision matrices that are estimated over the
#' same variables/features. This function visualizes in a single graph the
#' edges that are unique to the respective classes. Hence, it gives the
#' differential graph. Edges unique to \code{P1} are colored according to
#' \code{P1color}. Edges unique to \code{P2} are colored according to
#' \code{P2color}. Dashed lines indicate negative precision elements while
#' solid lines indicate positive precision elements.
#' 
#' The default layout is according to the Fruchterman-Reingold algorithm
#' (1991). Most layout functions supported by \code{\link{igraph}} are
#' supported (the function is partly a wrapper around certain
#' \code{\link{igraph}} functions). The igraph layouts can be invoked by a
#' \code{character} that mimicks a call to a \code{\link{igraph}} layout
#' functions in the \code{lay} argument. When using \code{lay = NULL} one can
#' specify the placement of vertices with the \code{coords} argument. The row
#' dimension of this matrix should equal the number of vertices. The column
#' dimension then should equal 2 (for 2D layouts) or 3 (for 3D layouts). The
#' \code{coords} argument can also be viewed as a convenience argument as it
#' enables one, e.g., to layout a graph according to the coordinates of a
#' previous call to \code{Ugraph}. If both the the lay and the coords arguments
#' are not \code{NULL}, the lay argument takes precedence.
#' 
#' @param P1 Sparsified precision \code{matrix} for class 1.
#' @param P2 Sparsified precision \code{matrix} for class 2.
#' @param lay A \code{character} mimicking a call to \code{\link{igraph}}
#' layout functions. Determines the placement of vertices.
#' @param coords A \code{matrix} containing coordinates. Alternative to the
#' lay-argument for determining the placement of vertices.
#' @param Vsize A \code{numeric} determining the vertex size.
#' @param Vcex A \code{numeric} determining the size of the vertex labels.
#' @param Vcolor A \code{character} (scalar or vector) determining the vertex
#' color.
#' @param VBcolor A \code{character} determining the color of the vertex
#' border.
#' @param VLcolor A \code{character} determining the color of the vertex
#' labels.
#' @param P1color A \code{character} determining the color of edges unique to
#' P1.
#' @param P2color A \code{character} determining the color of edges unique to
#' P2.
#' @param main A \code{character} giving the main figure title.
#' @return The function returns a graph.
#' @author Carel F.W. Peeters <cf.peeters@@vumc.nl>
#' @seealso \code{\link{Ugraph}}
#' @references Csardi, G. and Nepusz, T. (2006). The igraph software package
#' for complex network research. InterJournal, Complex Systems 1695.
#' http://igraph.sf.net
#' 
#' Fruchterman, T.M.J., and Reingold, E.M. (1991). Graph Drawing by
#' Force-Directed Placement. Software: Practice & Experience, 21: 1129-1164.
#' @examples
#' 
#' ## Obtain some (high-dimensional) data, class 1
#' p = 25
#' n = 10
#' set.seed(333)
#' X = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X)[1:25] = letters[1:25]
#' 
#' ## Obtain some (high-dimensional) data, class 2
#' set.seed(123456)
#' X2 = matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X2)[1:25] = letters[1:25]
#' 
#' ## Obtain regularized precision under optimal penalty, classes 1 and 2
#' OPT  <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)
#' OPT2 <- optPenalty.LOOCV(X2, lambdaMin = .5, lambdaMax = 30, step = 100)
#' 
#' ## Determine support regularized standardized precision under optimal penalty
#' PC0  <- sparsify(symm(OPT$optPrec), threshold = "localFDR")$sparseParCor
#' PC02 <- sparsify(symm(OPT2$optPrec), threshold = "localFDR")$sparseParCor
#' 
#' ## Visualize differential graph
#' DiffGraph(PC0, PC02)
#' 
#' @export DiffGraph
DiffGraph <- function(P1, P2, lay = "layout_with_fr", coords = NULL,
                      Vsize = 15, Vcex = 1, Vcolor = "orangered",
                      VBcolor = "darkred", VLcolor = "black",
                      P1color = "red", P2color = "green", main = ""){
  ##############################################################################
  # - Function that visualized a differential graph, i.e., the edges that are
  #   unique for 2 class-specific graphs over the same vertices
  # - Function is partly a wrapper around certain 'igraph' functions
  # - P1      > Sparsified precision matrix for class 1
  # - P2      > Sparsified precision matrix for class 2
  # - lay     > determines layout of the graph. Most layouts in 'layout{igraph}'
  #             are accepted. Default = layout_with_fr.
  # - coords  > matrix of coordinates to determine layout of the graph.
  #             The row dimension should equal the number of (pruned) vertices.
  #             The column dimension should equal 2 (for 2D layouts) or
  #             3 (for 3D layouts). Enables one, e.g., to layout the graph
  #             according to the coordinates of a previous call to Ugraph.
  #             If both the the lay and the coords arguments are not NULL,
  #             the lay argument takes precedence
  # - Vsize   > gives vertex size, default = 15
  # - Vcex    > gives size vertex labels, default = 1
  # - Vcolor  > gives vertex color, default = "orangered", must be character.
  #             May also be a character vector.
  # - VBcolor > gives color of the vertex border, default = "darkred"
  # - VLcolor > gives color of the vertex labels, default = "black"
  # - P1color > gives color of edges unique to P1, default = "red"
  # - P2color > gives color of edges unique to P1, default = "green"
  # - main    > character specifying heading figure, default = ""
  #
  # NOTES
  # - Will be expanded with legend options
  ##############################################################################

  # Dependencies
  # require("igraph")
  # require("reshape")

  if (!is.matrix(P1)){
    stop("Input (P1) should be a matrix")
  }
  else if (nrow(P1) != ncol(P1)){
    stop("Input (P1) should be square matrix")
  }
  else if (!is.matrix(P2)){
    stop("Input (P2) should be a matrix")
  }
  else if (nrow(P2) != ncol(P2)){
    stop("Input (P2) should be square matrix")
  }
  else if (nrow(P1) != nrow(P2)){
    stop("Inputs (P1 and P2) should be of the same dimensions")
  }
  else if (!((length(intersect(lay,c("layout_as_star", "layout_as_tree",
                                     "layout_in_circle", "layout_nicely",
                                     "layout_with_dh", "layout_with_gem",
                                     "layout_with_graphopt", "layout_on_grid",
                                     "layout_with_mds", "layout_components",
                                     "layout_on_sphere", "layout_randomly",
                                     "layout_with_fr", "layout_with_kk",
                                     "layout_with_lgl"))) > 0) | is.null(lay))){
    stop("lay should be 'NULL' or one of
         {'layout_as_star', 'layout_as_tree',
         'layout_in_circle', 'layout_nicely',
         'layout_with_dh', 'layout_with_gem',
         'layout_with_graphopt', 'layout_on_grid',
         'layout_with_mds', 'layout_components',
         'layout_on_sphere', 'layout_randomly',
         'layout_with_fr', 'layout_with_kk',
         'layout_with_lgl'}")
  }
  else if (!is.null(coords) & !inherits(coords, "matrix")){
    stop("Input (coords) is of wrong class")
  }
  else if (is.null(lay) & is.null(coords)){
    stop("Input (lay) and input (coords) cannot be both NULL")
  }
  else if (class(Vsize) != "numeric"){
    stop("Input (Vsize) is of wrong class")
  }
  else if (length(Vsize) != 1){
    stop("Length Vsize must be one")
  }
  else if (Vsize <= 0){
    stop("Vsize must be positive")
  }
  else if (class(Vcex) != "numeric"){
    stop("Input (Vcex) is of wrong class")
  }
  else if (length(Vcex) != 1){
    stop("Length Vcex must be one")
  }
  else if (Vcex <= 0){
    stop("Vcex must be positive")
  }
  else if (class(Vcolor) != "character"){
    stop("Input (Vcolor) is of wrong class")
  }
  else if (length(Vcolor) != 1 & length(Vcolor) != nrow(P1)){
    stop("Length Vcolor must be either one
         or equal to row (or column) dimension of P1 (P2)")
  }
  else if (class(VBcolor) != "character"){
    stop("Input (VBcolor) is of wrong class")
  }
  else if (length(VBcolor) != 1){
    stop("Length VBcolor must be one")
  }
  else if (class(VLcolor) != "character"){
    stop("Input (VLcolor) is of wrong class")
  }
  else if (length(VLcolor) != 1){
    stop("Length VLcolor must be one")
  }
  else if (class(P1color) != "character"){
    stop("Input (P1color) is of wrong class")
  }
  else if (length(P1color) != 1){
    stop("Length P1color must be one")
  }
  else if (class(P2color) != "character"){
    stop("Input (P2color) is of wrong class")
  }
  else if (length(P2color) != 1){
    stop("Length P2color must be one")
  }
  else if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  else {
    # Preliminaries 1
    AM1     <- adjacentMat(P1)
    AM2     <- adjacentMat(P2)
    DiffMat <- AM1 - AM2

    # Preliminaries 2
    WeightMat <- DiffMat
    WeightMat[WeightMat == 1]  <- P1[WeightMat == 1]
    WeightMat[WeightMat == -1] <- P2[WeightMat == -1]

    # Preliminaries 3
    AM <- adjacentMat(DiffMat)
    GA <- graph.adjacency(AM, mode = "undirected")

    # Layout specification
    if(is.null(lay)){
      if(dim(coords)[1] != length(V(GA))){
        stop("Row dimension of input (coords) does not match the
             number of vertices to be plotted")
      } else if (dim(coords)[2] > 3){
        stop("Column dimension of input (coords) exceeds the number
             of dimensions that can be visualized")
      } else {lays = coords}
      }
    else{
      if(lay == "layout_as_star"){
        lays = igraph::layout_as_star(GA)}
      if(lay == "layout_as_tree"){
        lays = igraph::layout_as_tree(GA)}
      if(lay == "layout_in_circle"){
        lays = igraph::layout_in_circle(GA)}
      if(lay == "layout_nicely"){
        lays = igraph::layout_nicely(GA)}
      if(lay == "layout_with_dh"){
        lays = igraph::layout_with_dh(GA)}
      if(lay == "layout_with_gem"){
        lays = igraph::layout_with_gem(GA)}
      if(lay == "layout_with_graphopt"){
        lays = igraph::layout_with_graphopt(GA)}
      if(lay == "layout_on_grid"){
        lays = igraph::layout_on_grid(GA)}
      if(lay == "layout_with_mds"){
        lays = igraph::layout_with_mds(GA)}
      if(lay == "layout_components"){
        lays = igraph::layout_components(GA)}
      if(lay == "layout_on_sphere"){
        lays = igraph::layout_on_sphere(GA)}
      if(lay == "layout_randomly"){
        lays = igraph::layout_randomly(GA)}
      if(lay == "layout_with_fr"){
        lays = igraph::layout_with_fr(GA)}
      if(lay == "layout_with_kk"){
        lays = igraph::layout_with_kk(GA)}
      if(lay == "layout_with_lgl"){
        lays = igraph::layout_with_lgl(GA)}
    }

    # Visualize
    Names <- colnames(DiffMat)
    colnames(DiffMat) = rownames(DiffMat) <- seq(1, ncol(DiffMat), by = 1)
    colnames(WeightMat) = rownames(WeightMat) <- seq(1, ncol(WeightMat), by = 1)
    Mmelt <- melt(DiffMat)
    Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
    Mmelt <- Mmelt[Mmelt$value != 0,]
    Mmelt2 <- melt(WeightMat)
    Mmelt2 <- Mmelt2[Mmelt2$X1 > Mmelt2$X2,]
    Mmelt2 <- Mmelt2[Mmelt2$value != 0,]
    E(GA)$weight <- Mmelt$value
    E(GA)$width  <- abs(E(GA)$weight)*1.5
    E(GA)$color  <- P1color
    E(GA)[E(GA)$weight == -1]$color <- P2color
    E(GA)$weight <- Mmelt2$value
    E(GA)[E(GA)$weight < 0]$style <- "dashed"
    E(GA)[E(GA)$weight > 0]$style <- "solid"
    plot(GA,
         layout = lays,
         vertex.size = Vsize,
         vertex.label.family = "sans",
         vertex.label.cex = Vcex,
         vertex.color = Vcolor,
         vertex.frame.color = VBcolor,
         vertex.label.color = VLcolor,
         edge.color = E(GA)$color,
         edge.lty = E(GA)$style,
         main = main)
  }
}




################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module D: rags2ridges DCMG
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# In development






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module E: rags2ridges Robust
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# In development






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Deprecated Collection
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# See R/rags2ridgesDepr.R
