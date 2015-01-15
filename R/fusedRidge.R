
.fusedUpdate <- function(k0, OmegaList, SList, targetList, ns, lambda1, lambda2) {
  K <- length(SList)
  stopifnot(k0 %in% seq_len(K))

  # Modified S
  b <- lambda2/ns[k0]
  tmp <- lapply(setdiff(seq_len(K), k0),
                function(k) OmegaList[[k]] - targetList[[k]])
  SS <- SList[[k0]] - b*Reduce('+', tmp)

  # Modified lambda
  a <- (lambda1 + 2*(K-1)*lambda2)/ns[k0]
  return(ridgeS(SS, lambda = a, target = targetList[[k0]]))
}

fusedRidgeS <- function(SList, targetList, ns, lambda1, lambda2,
                        max.ite = 100, verbose = TRUE, eps = 1e-6) {
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
