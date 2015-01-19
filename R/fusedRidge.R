

createS <- function(n, p) {
  ##############################################################################
  # - Simulate some random symmetric square matrices from uncorrolated noise
  # - n > A vector of sample sizes
  # - p > An integer giving the dimension
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
  if (K == 1) ans <- ans[[1]]  # Simplify output is ns is length 1
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


.fusedUpdate <- function(k0, PList, SList, TList, ns, lambda1, lambda2) {
  ##############################################################################
  # - (Internal) Update the scatter matrices
  # - k0      > An integer giving the class estimate to be updated.
  # - PList   > A list of matrices giving the current precision estimates.
  # - SList   > A list of sample correlation matrices the same size as PList.
  # - TList   > A list of target matrices the same size as SList
  # - ns      > A vector giving the sample sizes.
  # - lambda1 > The ridge penalty (a postive number).
  # - labmda2 > The fused penalty (a non-negative number).
  ##############################################################################

  b <- lambda2/ns[k0]                                          # Modded penalty2
  OmT <- mapply(`-`, PList[-k0], TList[-k0], SIMPLIFY = FALSE) # Omega - Target
  S0 <- SList[[k0]] - b*Reduce(`+`, OmT)                       # Modded S
  a <- 0.5*(lambda1 + 2*(length(SList) - 1)*lambda2)/ns[k0]    # Modded penalty1
  return(ridgeS(S0, lambda = a, target = TList[[k0]]))
}


fusedRidgeS <- function(SList, ns, TList = lapply(SList, default.target),
                        lambda1, lambda2,
                        max.ite = 100L, verbose = TRUE, eps = 1e-5) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda1 and lambda2
  # - SList   > A list of sample correlation matrices.
  # - TList   > A list as SList of target matrices with the same structure as
  #             SList. Default is given by default.target.
  # - ns      > A vector with the same length as SList giving the sample
  #             sizes.
  # - lambda1 > The ridge penalty (a postive number)
  # - lambda2 > The fused penalty (a non-negative number)
  # - max.ite > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - eps     > numeric. A positive convergence criterion.
  ##############################################################################

  stopifnot(length(SList) == length(TList))
  K <- length(SList)  # Number of groups
  PList <- SList      # Initialize estimates

  if (verbose) {
    cat("Iteration:  | Difference in Frobenious norm for k = ( 1 2 ... K )\n")
  }
  diffs <- rep(NA, K)
  for (i in seq_len(max.ite)) {
    for (k in seq_len(K)) {
      tmpOmega <- .fusedUpdate(k0 = k, PList = PList, SList = SList,
                               TList = TList, ns = ns,
                               lambda1 = lambda1, lambda2 = lambda2)

      if (!isSymmetric(tmpOmega)) {
        # Unfortunalty we need to do this check, the current isSymmetric is too
        # strict in the floating point precision.
        # (tol should be ~ sqrt(.Machine$double.eps) and not
        # 100*.Machine$double.eps)
        if (verbose) {cat("S")}
        tmpOmega <- symm(tmpOmega)
      } else {
        if (verbose) {cat(" ")}
      }

      diffs[k] <- .RelativeFrobeniusLoss(tmpOmega, PList[[k]])
      PList[[k]] <- tmpOmega
    }
    if (verbose) {
      cat(": i =", sprintf("%-2d", i), "| diffs = (", diffs, ")\n")
    }
    if (max(diffs) < eps) {
      break
    }
  }
  if (i == max.ite) {
    warning("Maximum iterations (", max.ite, ") hit")
  }
  return(PList)
}





optFusedPenalty.LOOCV <- function(YList,
                                  lambda1Min, lambda1Max, step1 = 100,
                                  lambda2Min = lambda1Min,
                                  lambda2Max = lambda1Max, step2 = step1,
                                  TList, verbose = TRUE) {
  if (missing(TList)) {  # If TList is not provided
    TList <- lapply(YList, function(Y) default.target(covML(Y)))
  }

  LLs     <- numeric()
  lambda1s <- exp(seq(log(lambda1Min), log(lambda1Max), length.out = step1))
  lambda2s <- exp(seq(log(lambda2Min), log(lambda2Max), length.out = step2))

  # Calculate CV scores
  if (verbose) {
    cat("Calculating cross-validated negative log-likelihoods...\n")
  }
  slh <- matrix(NA, step1, step2)
  total.n <- sum(sapply(YList, nrow))
  for (l1 in seq_along(lambda1s)) {
    for (l2 in seq_along(lambda2s)) {
      for (k in seq_along(YList)) {
        for (i in seq_len(nrow(YList[[k]]))) {

          SList <- lapply(YList, covML)
          SList[[k]] <- covML(YList[[k]][-i, , drop = FALSE])

          PList <- fusedRidgeS(SList, ns, TList,
                               lambda1 = lambda1s[l1],
                               lambda2 = lambda2s[l2], verbose = FALSE)

          SList[[k]] <- crossprod(YList[[k]][i, , drop = FALSE])
          slh[l1, l2] <- .FLL(SList, PList, ns)
          if (verbose){
            cat(sprintf("lambda1 = %.3f, lambda2 = %.3f, i = %-6d\n",
                        lambda1s[l1], lambda2s[l2], i, total.n))
          }
        }
      }
    }
  }
  return(list(lambda1s, lambda2s, slh))
}
