
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



.fusedUpdate <- function(k0, PList, SList, TList, ns, lambda1, LambdaP) {
  ##############################################################################
  # - (Internal) Update the scatter matrices
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



fusedRidgeS <- function(SList, ns, TList = lapply(SList, default.target),
                        lambda1, LambdaP, lambda2,
                        maxit = 100L, verbose = TRUE, eps = 1e-5) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda1 and lambda2
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
  # - eps     > numeric. A positive convergence criterion.
  ##############################################################################

  stopifnot(length(SList) == length(TList))
  K <- length(SList)  # Number of groups
  PList <- SList      # Initialize estimates

  if (!missing(LambdaP) && !missing(lambda2)) {
    stop("Supply only either LambdaP or lambda2.")
  } else if (missing(LambdaP) && missing(lambda2)) {
    stop("Either LambdaP or lambda2 must be given.")
  } else if (missing(LambdaP) && !missing(lambda2)) {
    LambdaP <- matrix(lambda2, K, K)
  }

  if (verbose) {
    cat("Iteration:  | Difference in Frobenious norm for k = ( 1 2 ... K )\n")
  }
  diffs <- rep(NA, K)
  for (i in seq_len(maxit)) {
    for (k in seq_len(K)) {
      tmpOmega <- .fusedUpdate(k0 = k, PList = PList, SList = SList,
                               TList = TList, ns = ns, lambda1 = lambda1,
                               LambdaP = LambdaP)

      diffs[k] <- .RelativeFrobeniusLoss(tmpOmega, PList[[k]])
      PList[[k]] <- tmpOmega
    }
    if (verbose) {
      cat("i =", sprintf("%-2d", i), "| diffs = (", diffs, ")\n")
    }
    if (max(diffs) < eps) {
      break
    }
  }
  if (i == maxit) {
    warning("Maximum iterations (", maxit, ") hit")
  }
  return(PList)
}


.fcvl <- function(lambdas, YList, TList) {
  ##############################################################################
  # - Helper function to perform fused LOOCV
  # - lambdas > A vector of penalties where lambdas[1] is lambda1 and
  #             lambdas[2] is lambda2.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
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
                            lambda1 = lambdas[1], lambda2 = lambdas[2],
                            maxit = 1000, verbose = FALSE, eps = 1e-4)
      slh <- c(slh, .LL(Sik, PList[[k]]))
    }
  }
  return(mean(slh))
}



.afcvl <- function(lambdas, YList, TList) {
  ##############################################################################
  # - (Internal) For at given lambda1 and lambda2, compute the approximate
  # - LOOCV loos
  # - The complete penalty graph is assumed  here.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  ##############################################################################

  ns <- sapply(YList, nrow)
  K <- length(ns)
  SList <- lapply(YList, covML)
  PList <- fusedRidgeS(SList = SList, TList = TList, ns = ns,
                       lambda1 = lambdas[1], lambda2 = lambdas[2],
                       verbose = FALSE)
  n.tot <- sum(ns)
  nll <- rags2ridges:::.FLL(SList = SList, PList = PList, ns)/n.tot
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



optFusedPenalty.LOOCV <- function(YList,
                                  lambda1Min, lambda1Max,
                                  step1 = 20,
                                  lambda2Min = lambda1Min,
                                  lambda2Max = lambda1Max,
                                  step2 = step1,
                                  TList,
                                  approximate = FALSE,
                                  verbose = TRUE) {
  ##############################################################################
  # - Simple leave one-out cross validation for the fused ridge estimator
  # - on a grid to determine optimal lambda1 and lambda2.
  # - The complete penalty graph is used here.
  # - YList      > A list of length K of matrices of observations with samples
  #                in the rows and variables in the columns.
  # - lambda1Min > Start lambda1 value, the ridge penalty
  # - lambda1Max > End lambda1 value
  # - step1      > Number of evaluations
  # - lambda2Min > As lambda1Min for the fused penalty. Default is lambda1Min.
  # - lambda2Max > As lambda1Max for the fused penalty. Default is lambda1Max.
  # - step2      > As step1 for the fused penalty. Default is step1.
  # - TList      > A list of length K of target matrices the same size
  #                as those of PList. Default is given by default.target.
  # - approximate > Should approximate LOOCV be used? Defaults to FALSE.
  #                 Approximate LOOCV is much faster.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
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
        cvl(lambdas = c(lambda1s[l1], lambda2s[l2]),
            YList = YList, TList = TList)
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
                                      approximate = FALSE,
                                      verbose = TRUE, ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. leave-one-out cross-validation
  # - using 2-dimensional BFGS optimization.
  # - The complete penalty graph is used here.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - approximate > logical. Should approximate LOOCV be used?
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - ...     > arguments passed to optim.
  ##############################################################################

  if (approximate) {  # Determine what function to use
    cvl <- .afcvl
  } else {
    cvl <- .fcvl
  }

  efcvl <- function(x) {  # Local reparameterized function
    cvl(exp(x), tYList, tTList)
  }

  # Get sensible starting value for lambda1 (choosing lambda2 to be zero)
  st <- optimize(function(x) efcvl(c(x, 0)), lower = 0, upper = 50)

  # Start at lambda2 point 0
  ans <- optim(c(st$minimum, 0), fn = efcvl, ...,
               method = "BFGS", control = list(trace = verbose))

  return(ans)
}





