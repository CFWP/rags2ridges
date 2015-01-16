

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


.symmetrize <- function(X) {
  ##############################################################################
  # - Ensure a matrix is symmetric by the average
  # - X > A square matrix
  ##############################################################################

  return((X + t(X))/2)
}


.FLL <- function(SList, PList, n){
  ##############################################################################
  # - Function that computes the value of the (negative) combined log-likelihood
  # - Slist > A list sample covariance matrices for each class
  # - Plist > A list of the same length as (Slist) of precision matrices
  #  (possibly regularized inverse of covariance or correlation matrices)
  # - n > A vector of sample sizes of the same length as Slist.
  ##############################################################################

  LLs <- mapply(.LL, SList, PList)
  return(sum(n*LLs))
}


.fusedUpdate <- function(k0, OList, SList, TList,
                         ns, lambda1, lambda2) {
  ##############################################################################
  # - (Internal) Update the scatter matrices
  # - k0      > An integer giving the class estimate to be updated.
  # - OList   > A list of matrices giving the current estimate (Omega).
  # - SList   > A list of sample correlation matrices the same size as OList.
  # - TList   > A list of target matrices the same size as SList
  # - ns      > A vector giving the sample sizes.
  # - lambda1 > The ridge penalty (a postive number).
  # - labmda2 > The fused penalty (a non-negative number).
  ##############################################################################

  # Modify sample covariance matrices
  b <- lambda2/ns[k0]                                    # Modified penalty2
  OmT <- mapply(`-`, OList[-k0], TList[-k0])    # Omega minus Target
  S0 <- SList[[k0]] - b*Reduce(`+`, OmT)                 # Modified S

  # Modified lambda
  a <- (lambda1 + 2*(length(SList) - 1)*lambda2)/ns[k0]  # Modified penalty2
  return(ridgeS(S0, lambda = a, target = TList[[k0]]))
}


fusedRidgeS <- function(SList, TList, ns, lambda1, lambda2,
                        max.ite = 100L, verbose = TRUE, eps = 1e-7) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda1 and lambda2
  # - SList   > A list of sample correlation matrices.
  # - TList   > A list as SList of target matrices with the same structure as
  #             SList.
  # - ns      > A vector with the same length as SList giving the sample
  #             sizes.
  # - lambda1 > The ridge penalty (a postive number)
  # - labmda2 > The fused penalty (a non-negative number)
  # - max.ite > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - eps     > numeric. A positive convergence criterion.
  ##############################################################################

  K <- length(SList)  # Number of groups
  OList <- SList  # Initialize estimates

  if (verbose) {
    cat("Iteration:  | Diff. in Omega in Frobenious norm for k = 1, ..., K\n")
  }
  diffs <- rep(NA, K)
  for (i in seq_len(max.ite)) {
    for (k in seq_len(K)) {
      tmpOmega <- .fusedUpdate(k0 = k, OList = OList, SList = SList,
                               TList = TList, ns = ns,
                               lambda1 = lambda1, lambda2 = lambda2)

      if (!isSymmetric(tmpOmega)) {
        # Unfortunalty we need to do this check, the current isSymmetric is too
        # strict in the floating point precision.
        # (tol should be ~ sqrt(.Machine$double.eps) and not
        # 100*.Machine$double.eps)
        if (verbose) {cat("S")}
        tmpOmega <- .symmetrize(tmpOmega)
      } else {
        if (verbose) {cat(" ")}
      }

      diffs[k] <- .FrobeniusLoss(tmpOmega, OList[[k]])/sum(tmpOmega^2)
      OList[[k]] <- tmpOmega
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
  return(OList)
}





