
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


symmetrize <- function(X) {
  ##############################################################################
  # - Ensure a matrix is symmetric by the average
  # - X > A square matrix
  ##############################################################################

  return((X + t(X))/2)
}


.fusedUpdate <- function(k0, OmegaList, SList, targetList,
                         ns, lambda1, lambda2) {
  ##############################################################################
  # - (Internal) Update the scatter matrices
  # - k0         > An integer giving the class estimate to be updated.
  # - OmegaList  > A list of matrices giving the current estimate.
  # - SList      > A list of sample correlation matrices the
  #                same size as OmegaList.
  # - targetList > A list of target matrices the same size as SList
  # - ns         > A vector giving the sample sizes.
  # - lambda1    > The ridge penalty (a postive number).
  # - labmda2    > The fused penalty (a non-negative number).
  ##############################################################################

  stopifnot(k0 %in% seq_len(K))

  # Modify sample covariance matrices
  b <- lambda2/ns[k0]
  OmT <- mapply(`-`, OmegaList[-k0], targetList[-k0])
  SS <- SList[[k0]] - b*Reduce(`+`, OmT)

  # Modified lambda`
  a <- (lambda1 + 2*(K-1)*lambda2)/ns[k0]
  return(ridgeS(SS, lambda = a, target = targetList[[k0]]))
}


fusedRidgeS <- function(SList, targetList, ns, lambda1, lambda2,
                        max.ite = 100, verbose = TRUE, eps = 1e-6) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda1 and lambda2
  # - SList      > A list of sample correlation matrices.
  # - targetList > A list as SList of target matrices with
  #                the same structure as SList.
  # - ns         > A vector with the same length as SList giving the sample
  #                sizes.
  # - lambda1    > The ridge penalty (a postive number)
  # - labmda2    > The fused penalty (a non-negative number)
  ##############################################################################

  K <- length(SList)

  # Initialize (corresponds to lambda2 = 0)
  #OmegaList <- lapply(SList, ridgeS, lambda = lambda1)
  #OmegaList <- lapply(seq_len(K), function(i) ridgeS(SList[[i]],lambda1/ns[i]))
  OmegaList <- SList

  if (verbose) {
    cat("Iteration:  | Diff. in Omega in Frobenious norm for k = 1, ..., K\n")
  }
  diffs <- rep(NA, K)
  for (i in seq_len(max.ite)) {
    for (k in seq_len(K)) {
      tmpOmega <- .fusedUpdate(k = k, OmegaList = OmegaList, SList = SList,
                               targetList = targetList, ns = ns,
                               lambda1 = lambda1, lambda2 = lambda2)

      if (!isSymmetric(tmpOmega)) { # Unfortunalty we need to do this
        if (verbose) cat("S")
        tmpOmega <- symmetrize(tmpOmega)
      } else {
        if (verbose) cat(" ")
      }

      diffs[k] <- .FrobeniusLoss(tmpOmega, OmegaList[[k]])
      OmegaList[[k]] <- tmpOmega
    }
    if (verbose) cat(": i =", sprintf("%-2d", i), "| diffs = (", diffs, ")\n")
    if (max(diffs) < eps) {
      break
    }
  }
  if (i == max.ite) {
    warning("Maximum iterations (", max.ite, ") hit")
  }
  return(OmegaList)
}
