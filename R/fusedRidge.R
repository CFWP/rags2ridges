
createS <- function(n, p) {
  ##############################################################################
  # - Simulate some random symmetric square matrices from uncorrelated noise
  # - n > A vector of sample sizes
  # - p > An integer giving the dimension
  # - Returns a list of matrices if n has length greater than 1
  ##############################################################################

  K <- length(n)

  # Construct names
  c <- which.max(p <= 25^(1:5))
  x <- expand.grid(rep(list(LETTERS), c))
  nms <- do.call(paste0, x)

  # Construct list
  ans <- list()
  for (i in seq_len(K)) {
    ans[[i]] <- covML(matrix(rnorm(n[i]*p), nrow = n[i], ncol = p))
    if (p <= 17576) {
      colnames(ans[[i]]) <- rownames(ans[[i]]) <- nms[1:p]
    }
  }
  if (K == 1) {ans <- ans[[1]]}  # Simplify output is ns is length 1
  return(ans)
}



.FLL <- function(SList, PList, ns){
  ##############################################################################
  # - Function that computes the value of the (negative) combined log-likelihood
  # - Slist > A list sample covariance matrices for each class
  # - Plist > A list of the same length as (Slist) of precision matrices
  #   (possibly regularized inverse of covariance or correlation matrices)
  # - ns > A vector of sample sizes of the same length as Slist.
  ##############################################################################

  LLs <- mapply(.LL, SList, PList)
  return(sum(ns*LLs))
}

.PFLL <- function(SList, PList, ns, TList, lambda1, LambdaP){
  ##############################################################################
  # - Function that computes the value of the (negative) penalized combined
  #   log-likelihood
  # - Slist   > A list sample covariance matrices for each class
  # - Plist   > A list of the same length as (Slist) of precision matrices
  #             (possibly regularized inverse of covariance  matrices)
  # - ns      > A vector of sample sizes of the same length as Slist.
  # - TList   > List of target matrices
  # - lambda1 > ridge penalty
  # - LambdaP > fused penalty matrix
  ##############################################################################

  penalty <- 0
  for (k1 in seq_along(SList)) {
    for (k2 in seq_len(k1)) {
      if (k1 == k2) { # Ridge penalty
        penalty <- penalty +
          lambda1*.FrobeniusLoss(SList[[k1]], TList[[k1]])
      } else {  # Fused contribution
        penalty <- penalty +
          LambdaP[k1, k2]*.FrobeniusLoss(SList[[k1]] - TList[[k1]],
                                         SList[[k2]] - TList[[k2]])
      }
    }
  }
  penalty <- penalty/2

  ans <- .FLL(SList, PList, ns) + penalty
  return(ans)
}



.fusedUpdate <- function(k0, PList, SList, TList, ns, lambda1, LambdaP) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate.
  # - k0      > An integer giving the class estimate to be updated.
  # - PList   > A list of length K of matrices giving the current precision
  #             estimates.
  # - SList   > A list of length K of sample correlation matrices the same size
  #             as those of PList.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList
  # - ns      > A vector of length K giving the sample sizes.
  # - lambda1 > The ridge penalty (a postive number).
  # - LambdaP > A K by K symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where LambdaP[k1, k2] determine
  #             the (rate of) shrinkage between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  ##############################################################################

  diag(LambdaP) <- 0  # Make sure there's zeros in the diagonal
  a <- (sum(LambdaP[k0, ]) + lambda1)/ns[k0]

  OmT <- mapply(`-`, PList[-k0], TList[-k0], SIMPLIFY = FALSE) # Omega - Target
  OmT <- mapply(`*`, LambdaP[k0, -k0]/ns[k0], OmT, SIMPLIFY = FALSE)
  S0 <- SList[[k0]] - Reduce(`+`, OmT)
  return(ridgeSArma(S0, lambda = a, target = TList[[k0]]))
}



.fusedUpdate2 <- function(k0, PList, SList, TList, ns, lambda1, LambdaP) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate.
  # - k0      > An integer giving the class estimate to be updated.
  # - PList   > A list of length K of matrices giving the current precision
  #             estimates.
  # - SList   > A list of length K of sample correlation matrices the same size
  #             as those of PList.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList
  # - ns      > A vector of length K giving the sample sizes.
  # - lambda1 > The ridge penalty (a postive number).
  # - LambdaP > A K by K symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where LambdaP[k1, k2] determine
  #             the (rate of) shrinkage between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  ##############################################################################

  lambdaa <- (lambda1 + sum(LambdaP[k0, -k0]))/ns[k0]

  p <- nrow(PList[[1]])
  M <- matrix(0, p, p)
  for (k in setdiff(seq_along(PList), k0)) {
    M <- M + (LambdaP[k, k0]/ns[k0])*(PList[[k]] - TList[[k]])
  }
  stopifnot(isSymmetric(M))
  return(ridgeSArma(SList[[k0]] - M, lambda = lambdaa, target = TList[[k0]]))
}



fusedRidgeS <- function(SList, ns, TList = lapply(SList, default.target),
                        lambda1, LambdaP, lambda2, PList,
                        maxit = 100L, verbose = TRUE, eps = 1e-4) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda1 and LambdaP
  # - SList   > A list of length K of sample correlation matrices the same size
  #             as those of PList.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ns      > A vector of length K giving the sample sizes.
  # - lambda1 > The ridge penalty (a postive number)
  # - LambdaP > The K by K symmetric adjacency matrix fused penalty graph
  #             with non-negative entries where LambdaP[k1, k2] determine the
  #             retainment of similarities between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  # - lambda2 > The non-negative fused penalty. Alternative to LambdaP if
  #             all pairwise penalties equal.
  # - maxit   > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - eps     > numeric. A positive convergence criterion. Default is 1e-4.
  ##############################################################################

  stopifnot(length(SList) == length(TList))
  K <- length(SList)  # Number of groups

  # Initialize estimates with the regular ridge for each class
  if (missing(PList)) {
    Spool <- Reduce(`+`, mapply("*", ns, SList, SIMPLIFY = FALSE))/sum(ns)
    PList <- list()
    for (i in seq_len(K)) {
      PList[[i]] <- ridgeSArma(Spool, lambda = lambda1, target = TList[[i]])
    }
  }
  stopifnot(length(SList) == length(PList))

  # Initialize estimates with the regular ridge for each class
  PList <- mapply(ridgeS, SList, lambda1, SIMPLIFY = FALSE)

  if (!missing(LambdaP) && !missing(lambda2)) {
    stop("Supply only either LambdaP or lambda2.")
  } else if (missing(LambdaP) && missing(lambda2)) {
    stop("Either LambdaP or lambda2 must be given.")
  } else if (missing(LambdaP) && !missing(lambda2)) {
    LambdaP <- matrix(lambda2, K, K)
  }

  if (verbose) {
    cat("Iter:   | difference in Frobenius norm        | -penalized log-lik\n")
    cat("init    | diffs = (", sprintf("%11e", rep(NA, K)), ")")
    cat(sprintf(" | -pll = %g\n", .PFLL(SList,PList,ns,TList,lambda1,LambdaP)))
  }


  tmpPList <- list()
  diffs <- rep(NA, K)
  i <- 1
  while (i <= maxit) {
    for (k in seq_len(K)) {
      tmpPList[[k]] <- .fusedUpdate(k0 = k, PList = PList, SList = SList,
                                    TList = TList, ns = ns, lambda1 = lambda1,
                                    LambdaP = LambdaP)
      diffs[k] <- .FrobeniusLoss(tmpPList[[k]], PList[[k]])
      PList[[k]] <- tmpPList[[k]]
    }

    if (verbose) {
      cat(sprintf("i = %-3d", i), "| diffs = (", sprintf("%6.5e", diffs), ")")
      cat(sprintf(" | -pll = %g\n",.PFLL(SList,PList,ns,TList,lambda1,LambdaP)))
    }

    if (max(diffs) < eps) {
      break
    }

    i <- i + 1
  }

  if (i == maxit + 1) {
    warning("Maximum iterations (", maxit, ") hit")
  }
  return(PList)
}



################################################################################
## -------------------------------------------------------------------------- ##
##
## Leave-one-out cross-validation (and approximation) in the fused setting
##
## -------------------------------------------------------------------------- ##
################################################################################



.fcvl <- function(lambda1, LambdaP, YList, TList, ...) {
  ##############################################################################
  # - A function to perform fused LOOCV
  # - lambda1 > numeric of length 1 giving ridge penalty.
  # - LambdaP > numeric matrix giving the fused penalty matrix.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ...     > Arguments passed to fusedRidgeS
  ##############################################################################

  K <- length(YList)
  ns.org <- sapply(YList, nrow)
  SList.org <- lapply(YList, covML)

  for (k in seq_len(K)) {
    nk <- nrow(YList[[k]])
    slh <- numeric()
    for (i in seq_len(nk)) {
      ns <- ns.org
      ns[k] <- ns[k] - 1
      SList <- SList.org
      SList[[k]] <- covML(YList[[k]][-i, , drop = FALSE])
      Sik    <- crossprod(YList[[k]][i,  , drop = FALSE])
      PList  <- fusedRidgeS(SList = SList, ns = ns, TList = TList,
                            lambda1 = lambda1, LambdaP = LambdaP,
                            verbose = FALSE, ...)
      slh <- c(slh, .LL(Sik, PList[[k]]))
    }
  }
  return(mean(slh))
}



.afcvl <- function(lambda1, LambdaP, YList, TList, ...) {
  ##############################################################################
  # - (Internal) For at given lambda1 and lambda2, compute the approximate
  # - LOOCV loss
  # - lambda1 > numeric of length 1 giving ridge penalty.
  # - LambdaP > numeric matrix giving the fused penalty matrix.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ...     > Arguments passed to fusedRidgeS
  ##############################################################################

  ns <- sapply(YList, nrow)
  K <- length(ns)
  SList <- lapply(YList, covML)
  PList <- fusedRidgeS(SList = SList, TList = TList, ns = ns,
                       lambda1 = lambda1, LambdaP = LambdaP,
                       verbose = FALSE, ...)
  n.tot <- sum(ns)
  nll <- .FLL(SList = SList, PList = PList, ns)/n.tot
  denom <- n.tot*(n.tot - 1)
  bias <- 0
  for (k in seq_along(ns)) {
    for (i in seq_len(ns[k])) {
      Sik <- crossprod(YList[[k]][i, , drop = FALSE])
      fac1 <- diag(nrow(Sik)) - Sik %*% PList[[k]]
      fac2 <- PList[[k]] %*% (SList[[k]] - Sik) %*% PList[[k]]
      bias <- bias  + sum(fac1 * fac2)/denom
    }
  }
  return(nll + bias)
}



.parseLambda <- function(Lambda) {
  ##############################################################################
  # - A function to parse a character matrix that defines the class of penalty
  #   graphs and unique parameters for cross validation. Returns a list of
  #   indices for each level to be penalized equally.
  #   This list is to be used to construct numeric matrices of penalties.
  # - Lambda > A square K by K character matrix defining the class penalty
  #            matrices to use. Entries with NA, "NA", "" (the empty string),
  #            or "0" are used specify that the pair should be omitted.
  ##############################################################################

  stopifnot(is.character(Lambda))
  stopifnot(is.matrix(Lambda))
  stopifnot(nrow(Lambda) == ncol(Lambda))

  Lambda[is.na(Lambda)] <- ""
  Lambda[Lambda %in% c("0", "NA")] <- ""

  lvls <- unique(c(Lambda))
  lvls <- lvls[lvls != ""]

  # For each non-empty level get boolean matrices
  ans <- lapply(lvls, function(lvl) which(lvl == Lambda, arr.ind = TRUE))
  names(ans) <- lvls
  return(ans)
}



.reconstructLambda <- function(lambdas, parsedLambda, K) {
  ##############################################################################
  # - Reconstruct the numeric penalty matrix Lambda from vector (lambdas) using
  #   output from .parseLambda output.
  # - lambdas      > A numeric vector of length K+1 where the first entry is the
  #                  ridge penalty and the remaining are the fused penalties.
  # - parsedLambda > A list of length K of matrix indicies.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  stopifnot(length(lambdas) == length(parsedLambda) + 1)
  LambdaP <- matrix(0, K, K)
  for (i in seq_along(parsedLambda)) {
    LambdaP[parsedLambda[[i]]] <- lambdas[i + 1]
  }
  return(LambdaP)
}



optFusedPenalty.LOOCV <- function(YList,
                                  lambda1Min, lambda1Max,
                                  step1 = 20,
                                  lambda2Min = lambda1Min,
                                  lambda2Max = lambda1Max,
                                  step2 = step1,
                                  TList,
                                  approximate = FALSE,
                                  ...,
                                  verbose = TRUE) {
  ##############################################################################
  #   Simple (approximate) leave one-out cross validation for the fused ridge
  #   estimator on a grid to determine optimal lambda1 and lambda2.
  #   The complete penalty graph is assumed.
  # - YList       > A list of length K of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - lambda1Min  > Start of lambda1 value, the ridge penalty
  # - lambda1Max  > End of lambda1 value
  # - step1       > Number of evaluations
  # - lambda2Min  > As lambda1Min for the fused penalty. Default is lambda1Min.
  # - lambda2Max  > As lambda1Max for the fused penalty. Default is lambda1Max.
  # - step2       > As step1 for the fused penalty. Default is step1.
  # - TList       > A list of length K of target matrices the same size
  #                 as those of PList. Default is given by default.target.
  # - approximate > Should approximate LOOCV be used? Defaults is FALSE.
  #                 Approximate LOOCV is much faster.
  # - ...         > Arguments passed to fusedRidgeS
  # - verbose     > logical. Print extra information. Defaults is TRUE.
  #
  # The function evaluates the loss on a log-equidistant grid.
  ##############################################################################

  if (missing(TList)) {  # If TList is not provided
    TList <- lapply(YList, function(Y) default.target(covML(Y)))
  }

  K <- length(YList)
  ns <- sapply(YList, nrow)

  # Choose lambdas log-equidistantly
  lambda1s <- exp(seq(log(lambda1Min), log(lambda1Max), length.out = step1))
  lambda2s <- exp(seq(log(lambda2Min), log(lambda2Max), length.out = step2))
  stopifnot(all(is.finite(lambda1s)))
  stopifnot(all(is.finite(lambda2s)))

  if (approximate) {
    cvl <- .afcvl
  } else {
    cvl <- .fcvl
  }

  # Calculate CV scores
  if (verbose) {
    cat("Calculating cross-validated negative log-likelihoods...\n")
  }

  slh <- matrix(NA, step1, step2)
  total.n <- sum(sapply(YList, nrow))
  for (l1 in seq_along(lambda1s)) {
    for (l2 in seq_along(lambda2s)) {
      slh[l1, l2] <-
        cvl(lambda1 = lambda1s[l1], LambdaP = matrix(lambda2s[l2], K, K),
            YList = YList, TList = TList, ...)
      if (verbose){
        cat(sprintf("lambda1 = %.3f (%d), lambda2 = %.3f (%d), -ll = %.3f\n",
                   lambda1s[l1],  l1, lambda2s[l2], l2, slh[l1, l2]))
      }
    }
  }
  return(list(lambda1 = lambda1s, lambda2 = lambda2s, fcvl = slh))
}



optFusedPenalty.LOOCVauto <- function(YList,
                                      TList,
                                      Lambda,
                                      approximate = FALSE,
                                      verbose = TRUE,
                                      maxit.fusedRidgeS = 1000,
                                      maxit.optim = 1000,
                                      optim.debug = FALSE,
                                      ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   leave-one-out cross-validation using multi-dimensional BFGS optimization
  #   routines.
  #
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - Lambda  > A K by K character matrix defining the class of penalty graph
  #             to use. The unique elements of Lambda specify the penalties to
  #             determine. Pairs can be left out using either of "", NA,
  #             "NA" or "0".
  # - approximate > logical. Should approximate LOOCV be used?
  # - verbose     > logical. Should the function print extra info. Defaults to
  #                 TRUE.
  # - maxit.fusedRidgeS > integer. Maximum number of iterations for fusedRidgeS
  # - maxit.optim       > integer. Maximum number of iterations for optim.
  # - optim.debug       > logical. If TRUE the raw output from optim is added
  #                       as an attribute to the output.
  # - ...               > arguments passed to optim.
  #
  # The function returns a list of length 4 with entries (1) lambda1,
  # (2) lambda2, (3) the optimal penalty matrix LambdaP, and (4) the value of
  # the loss in the optimum. If LambdaP is the complete graph, then lambda2 is
  # given. Otherwise lambda2 is NA.
  ##############################################################################

  K <- length(YList)

  # Handle Lambda
  if (missing(Lambda)) {
    Lambda <- matrix("A", K, K)
    diag(Lambda) <- ""
  }

  parsedLambda <- .parseLambda(Lambda)
  n.lambdas <- length(parsedLambda) + 1

  # Determine what loss function to use
  # lambdas[1] is the regular ridge penalty, while the remaning lambdas[-1]
  # correspond to the penalty matrix.
  # We also reparameterize to work on log-scale
  if (approximate) {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .afcvl(lambda1 = elambdas[1],
             LambdaP = .reconstructLambda(elambdas, parsedLambda, K),
             YList = YList, TList = TList, maxit = maxit.fusedRidgeS, ...)
    }
  } else {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .fcvl(lambda1 = elambdas[1],
            LambdaP = .reconstructLambda(elambdas, parsedLambda, K),
            YList = YList, TList = TList, maxit = maxit.fusedRidgeS, ...)
    }
  }

  # Get sensible starting value for lambda1 (choosing LambdaP to be zero)
  st <- optimize(function(x) cvl(c(x, rep(0, n.lambdas - 1))),
                 lower = -30, upper = 30)

  # Start LambdaP at 0
  lambdas.init <- c(st$minimum, rep(0, n.lambdas - 1))
  ans <- optim(lambdas.init, fn = cvl, ...,
               method = "BFGS",
               control = list(trace = verbose, maxit = maxit.optim))

  # Format optimal values
  opt.lambdas <- exp(ans$par)
  res <- list(lambda1 = opt.lambdas[1],
              lambda2 = NA,
              LambdaP = .reconstructLambda(opt.lambdas, parsedLambda, K),
              value = ans$value)
  lambda2 <- unique(opt.lambdas[-1])
  if (length(lambda2) == 1) {
    res$lambda2 <- lambda2
  }

  if (optim.debug) {
    attr(res, "optim.debug") <- ans
  }
  return(res)
}





################################################################################
## -------------------------------------------------------------------------- ##
##
## Automatic penalty matrix constructor
##
## -------------------------------------------------------------------------- ##
################################################################################



.charAdjMat <- function(fac, name = "X") {
  ##############################################################################
  # Create a character adjacency matrix from a factor
  # - fac  > A factor of some length.
  # - name > A character giving the text which should appear in the adjacent
  #          entries. If not a character, the object name is used.
  # Examples:
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "Y")
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = NULL)
  ##############################################################################

  p <- nlevels(fac)
  if (is.character(name)) {
    lab <- name
  } else {
    lab <- deparse(substitute(fac))
  }
  M <- matrix(lab, p, p)
  diag(M) <- ""
  dimnames(M) <- replicate(2, levels(fac), simplify = FALSE)
  return(M)
}

.char2num <- function(X) {
  ##############################################################################
  # Create a character adjacency matrix to a numeric one
  # - X  > A character matrix where "" signify non-adjacency.
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  # rags2ridges:::.char2num(A)
  ##############################################################################

  X[X != ""] <- "1"
  X[X == ""] <- "0"
  Y <- structure(as.numeric(X), dim = dim(X), dimnames = dimnames(X))
  return(Y)
}

.cartesianProd <- function(A, B) {
  ##############################################################################
  # Create a character adjacency matrix to a numeric one
  # - A     > A character matrix where "" signify non-adjacency.
  # - B     > A character matrix where "" signify non-adjacency.
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "X")
  # B <- rags2ridges:::.charAdjMat(factor(letters[4:5]), name = "Y")
  # rags2ridges:::.cartesianProd(A, B)
  ##############################################################################

  AI <- kronecker(.char2num(A), diag(nrow(B)), make.dimnames = TRUE)
  IB <- kronecker(diag(nrow(A)), .char2num(B), make.dimnames = TRUE)
  prod <- AI + IB

  ans <- kronecker(A, B, FUN = paste0, make.dimnames = TRUE)
  ans[!as.logical(prod)] <- ""
  return(ans)
}


default.penalty <- function(K, df,
                            type = c("Complete", "CartesianEqual",
                                     "CartesianUnequal")) {
  ##############################################################################
  # Select a one of standard penalty matrix types from a dataframe
  # - K     > The number of classes. Can also be list of length K such as
  #           the usual argument "SList".
  #           Can be omitted if 'df' is given.
  # - df    > A data.frame with K rows and with factors in the columns.
  #           Columns of type character are coerced to factors.
  #           Can be omitted when 'type == "Complete"'.
  # - type  > A character giving the type of fused penalty graph to construct.
  #           Should be 'Complete' (default), 'CartesianEqual', or
  #           'CartesianUnequal'.
  # Setting type == 'Complete' is the complete penalty graph with equal
  # penalties.
  # Setting type == 'CartesianEqual' corresponds to a penalizing along each
  # "direction" of factors with a common penalty.
  # Setting type == 'CartesianUnequal' corresponds to a penalizing each
  # direction of factors with individual penalties.
  ##############################################################################
  type <- match.arg(type)

  if (is.data.frame(K)) {
    df <- K
    K <- nrow(df)
  }

  if (missing(K) && !missing(df)) {
    K <- nrow(df)
  }

  if (is.list(K)) {
    K <- length(K)
  }

  if (missing(df)) {
    if (type != "Complete") {
      warning("No data.frame 'df' given and 'type' does not equal 'Complete'.",
              " Setting 'type' to 'Complete'")
      type <- "Complete"
    }
    df <- data.frame(Class = factor(seq_len(K)))
  }

  if (!all(sapply(df, is.factor))) {
    stop("Not all columns in the data.frame 'df' are factors")
  }

  stopifnot(K == nrow(df))

  if (type == "Complete") {

    M <- matrix("f", K, K)
    rownames(M) <- colnames(M) <- Reduce(":", df)
    diag(M) <- ""
    return(M)

  } else if (type == "CartesianEqual" || type == "CartesianUnequal") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.cartesianProd, adj.mats)

    if (type == "CartesianEqual") {
      M[M != ""] <- "f"
    }
    return(M)

  } else {
    stop("type", type, "not implemented yet!")
  }
}




