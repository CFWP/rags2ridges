
createS <- function(n, p,
                    topology = "identity",
                    dataset = FALSE,
                    precision = FALSE,
                    nonzero = 0.25,
                    m = 1L,
                    banded.n = 2L) {
  ##############################################################################
  # - Simulate some random symmetric square matrices from uncorrelated noise
  #   or datasets
  # - n          > A vector of sample sizes
  # - p          > An integer giving the dimension. p should be greater than 2.
  # - topology   > character. Specify the topology to simulate data from.
  #                See details.
  # - dataset    > logical. Should dataset instead of its sample covariance be
  #                returned? Default is FALSE.
  # - precision  > logical. Should the constructed precision matrix be returned?
  # - nonzero    > numeric of length 1 giving the value of the non-zero entries
  #                for some topologies
  # - m          > integer. The number of conditionally independent subgraphs.
  #                I.e. the number of blocks.
  # - banded.n   > interger. The number of off-diagonal bands used if
  #                topology is "banded". Use as paremeter if topology is
  #                "Watt-Strogatz".
  # - Returns a list of matrices if length(n) > 1.  The output is simplified if
  #   n has length 1, and only the matrix is returned
  ##############################################################################

  stopifnot(p > 1)
  stopifnot(m >= 1)

  if (dataset && precision) {
    stop("dataset and precision cannot be TRUE at the same time.")
  }

  topology <- match.arg(topology,
                        c("identity", "star", "clique", "complete",
                          "chain", "banded", "Barabassi", "small-world",
                          "scale-free", "Watts-Strogatz", "random-graph",
                          "Erdos-Renyi"))
  K <- length(n)

  # Construct the precision matrix "constructor"
  if (topology == "identity") {

    submat <- function(p) diag(p)

  } else if (topology == "star") {

    submat <- function(p) {
      subP <- diag(p)
      subP[1, seq(2, p)] <- subP[seq(2, p), 1] <- 1/seq(2, p)
      return(subP)
    }

  } else if (topology == "chain") {

    submat <- function(p) {
      s <- seq_len(p - 1)
      subP <- diag(p)
      subP[cbind(s, s + 1)] <- subP[cbind(s + 1, s)] <- nonzero
      return(subP)
    }

  } else if (topology == "clique" || topology == "complete") {

    submat <- function(p) {
      subP <- diag(p)
      subP[lower.tri(subP)] <- subP[upper.tri(subP)]  <- nonzero
      return(subP)
    }

  } else if (topology == "banded") {

    submat <- function(p) {
      if (banded.n > p) {
        stop("The number of bands cannot exceed the dimension of each block")
      }
      subP <- diag(p)
      for (j in seq(1, banded.n)) {
        s <- seq_len(p - j)
        subP[cbind(s, s + j)] <- subP[cbind(s + j, s)] <- 1/(j + 1)
      }
      return(subP)
    }

  } else if (topology == "Barabassi" || topology == "scale-free") {

    submat <- function(p) {
      G <- barabasi.game(p, power = 1, directed = FALSE)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else if (topology == "Watts-Strogatz" || topology == "small-world") {

    submat <- function(p) {
      G <- watts.strogatz.game(1, p, banded.n, 0.05)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else if (topology == "Erdos-Renyi" || topology == "random-graph") {

    submat <- function(p) {
      G <- erdos.renyi.game(p, 1/p)
      adj <- get.adjacency(G, sparse = FALSE)
      return(diag(p) + nonzero*adj)
    }

  } else {

    stop("requested topology not implemented yet.")

  }

  # Construct the block split
  blocks <- split(seq_len(p), ceiling(m*seq_len(p)/p))

  # Fill in blocks to construct full precisions
  P <- diag(p)
  for (b in blocks) {
    P[b, b] <- submat(length(b))
  }
  if (rcond(P) < sqrt(.Machine$double.eps)) {  # Check condition number
    warning("The generated precision matrix has a very high condition number ",
            "and the generated data might be unreliable.")
  }
  S <- solve(P)

  # Construct names
  n.letters <- which.max(p <= 26^(1:3))
  x <- expand.grid(rep(list(LETTERS), n.letters))
  nms <- do.call(paste0, x)

  # Construct list to fill and iterate over all classes
  ans <- vector("list", K)
  names(ans) <- paste0("class", seq_len(K))
  for (i in seq_len(K)) {

    if (precision) {
      ans[[i]] <- P
    } else {
      ans[[i]] <- rmvnormal(n = n[i], mu = rep(0, p), sigma = S)
      if (!dataset) {
        ans[[i]] <- covML(ans[[i]])
      }
    }
    if (p <= 17576) {  # Only give names for "small" dimensions
      colnames(ans[[i]]) <- nms[1:p]
      if (!dataset) {
        rownames(ans[[i]]) <- nms[1:p]
      }
    }
  }

  if (K == 1) {  # Simplify output if ns is length 1
    ans <- ans[[1]]
  }

  return(ans)
}



pooledS <- function(SList, ns, mle = TRUE) {
  ##############################################################################
  # - Computes the pooled covariance estimate
  # - Slist > A list sample covariance matrices for each class
  # - ns    > A vector of sample sizes of the same length as Slist.
  # - mle   > logical. If TRUE the biased MLE is used. If FALSE, the biased
  #           corrected estimate is used.
  ##############################################################################
  ans <- armaPooledS(SList = SList, ns = ns, mle = as.numeric(mle))
  dimnames(ans) <- dimnames(SList[[1]])
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

.PFLL <- function(SList, PList, ns, TList, lambda, lambdaFmat){
  ##############################################################################
  # - Function that computes the value of the (negative) penalized combined
  #   log-likelihood
  # - Slist   > A list sample covariance matrices for each class
  # - Plist   > A list of the same length as (Slist) of precision matrices
  #             (possibly regularized inverse of covariance  matrices)
  # - ns      > A vector of sample sizes of the same length as Slist.
  # - TList   > List of target matrices
  # - lambda  > ridge penalty
  # - lambdaFmat > fused penalty matrix
  ##############################################################################

  penalty <- 0
  for (k1 in seq_along(SList)) {
    for (k2 in seq_len(k1)) {
      if (k1 == k2) { # Ridge penalty
        penalty <- penalty +
          lambda*.FrobeniusLoss(SList[[k1]], TList[[k1]])
      } else {  # Fused contribution
        penalty <- penalty +
          lambdaFmat[k1, k2]*.FrobeniusLoss(SList[[k1]] - TList[[k1]],
                                         SList[[k2]] - TList[[k2]])
      }
    }
  }
  penalty <- penalty/2

  ans <- .FLL(SList, PList, ns) + penalty
  return(ans)
}



.fusedUpdate <- function(k0, PList, SList, TList, ns, lambda, lambdaFmat) {
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
  # - lambda > The ridge penalty (a postive number).
  # - lambdaFmat > A K by K symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambdaFmat[k1, k2] determine
  #             the (rate of) shrinkage between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  ##############################################################################

  diag(lambdaFmat) <- 0  # Make sure there's zeros in the diagonal
  a <- (sum(lambdaFmat[k0, ]) + lambda)/ns[k0]

  OmT <- mapply(`-`, PList[-k0], TList[-k0], SIMPLIFY = FALSE) # Omega - Target
  OmT <- mapply(`*`, lambdaFmat[k0, -k0]/ns[k0], OmT, SIMPLIFY = FALSE)
  S0 <- SList[[k0]] - Reduce(`+`, OmT)
  return(armaRidgeS(S0, target = TList[[k0]], lambda = a))
}



.fusedUpdate2 <- function(k0, PList, SList, TList, ns, lambda, lambdaFmat) {
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
  # - lambda > The ridge penalty (a postive number).
  # - lambdaFmat > A K by K symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambdaFmat[k1, k2] determine
  #             the (rate of) shrinkage between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  ##############################################################################

  lambdaa <- (lambda + sum(lambdaFmat[k0, -k0]))/ns[k0]

  p <- nrow(PList[[1]])
  M <- matrix(0, p, p)
  for (k in setdiff(seq_along(PList), k0)) {
    M <- M + (lambdaFmat[k, k0]/ns[k0])*(PList[[k]] - TList[[k]])
  }
  stopifnot(isSymmetric(M))
  return(armaRidgeS(SList[[k0]] - M, target = TList[[k0]], lambda = lambdaa))
}



ridgeS.fused <- function(SList, ns, TList = default.target.fused(SList, ns),
                         lambda, lambdaFmat, lambdaF, PList,
                         maxit = 100L, verbose = TRUE, eps = 1e-4) {
  ##############################################################################
  # - The fused ridge estimate for a given lambda and lambdaFmat
  # - SList   > A list of length K of sample correlation matrices the same size
  #             as those of PList.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ns      > A vector of length K giving the sample sizes.
  # - lambda > The ridge penalty (a postive number)
  # - lambdaFmat > The K by K symmetric adjacency matrix fused penalty graph
  #             with non-negative entries where lambdaFmat[k1, k2] determine the
  #             retainment of similarities between estimates in classes
  #             corresponding to SList[k1] and SList[k1].
  # - lambdaF > The non-negative fused penalty. Alternative to lambdaFmat if
  #             all pairwise penalties equal.
  # - maxit   > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - eps     > numeric. A positive convergence criterion. Default is 1e-4.
  ##############################################################################

  stopifnot(length(SList) == length(TList))
  K <- length(SList)  # Number of groups

  # Initialize estimates with the regular ridges from the pooled covariance
  if (missing(PList)) {
    Spool <- Reduce(`+`, mapply("*", ns, SList, SIMPLIFY = FALSE))/sum(ns)
    PList <- list()
    for (i in seq_len(K)) {
      PList[[i]] <- armaRidgeS(Spool, target = TList[[i]], lambda = lambda)
    }
  }
  stopifnot(length(SList) == length(PList))

  if (!missing(lambdaFmat) && !missing(lambdaF)) {
    stop("Supply only either lambdaFmat or lambdaF.")
  } else if (missing(lambdaFmat) && missing(lambdaF)) {
    stop("Either lambdaFmat or lambdaF must be given.")
  } else if (missing(lambdaFmat) && !missing(lambdaF)) {
    lambdaFmat <- matrix(lambdaF, K, K)
  }

  if (verbose) {
    cat("Iter:   | difference in Frobenius norm        | -penalized log-lik\n")
    cat("init    | diffs = (", sprintf("%11e", rep(NA, K)), ")")
    cat(sprintf(" | -pll = %g\n",
                .PFLL(SList, PList, ns, TList, lambda, lambdaFmat)))
  }


  tmpPList <- list()
  diffs <- rep(NA, K)
  i <- 1
  while (i <= maxit) {
    for (k in seq_len(K)) {
      tmpPList[[k]] <- .fusedUpdate(k0 = k, PList = PList, SList = SList,
                                    TList = TList, ns = ns, lambda = lambda,
                                    lambdaFmat = lambdaFmat)
      diffs[k] <- .FrobeniusLoss(tmpPList[[k]], PList[[k]])
      PList[[k]] <- tmpPList[[k]]
    }

    if (verbose) {
      cat(sprintf("i = %-3d", i), "| diffs = (", sprintf("%6.5e", diffs), ")")
      cat(sprintf(" | -pll = %g\n",
                  .PFLL(SList,PList,ns,TList,lambda,lambdaFmat)))
    }

    if (max(diffs) < eps) {
      break
    }

    i <- i + 1
  }

  if (i == maxit + 1) {
    warning("Maximum iterations (", maxit, ") hit")
  }

  # Keep dimnames and names
  for (k in seq_along(SList)) {
    dimnames(PList[[k]]) <- dimnames(SList[[k]])
  }
  names(PList) <- names(SList)


  return(PList)
}



################################################################################
## -------------------------------------------------------------------------- ##
##
## Leave-one-out cross-validation (and approximation) in the fused setting
##
## -------------------------------------------------------------------------- ##
################################################################################



.fcvl <- function(lambda, lambdaFmat, YList, TList, ...) {
  ##############################################################################
  # - A function to perform fused LOOCV
  # - lambda > numeric of length 1 giving ridge penalty.
  # - lambdaFmat > numeric matrix giving the fused penalty matrix.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ...     > Arguments passed to ridgeS.fused
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
      PList  <- ridgeS.fused(SList = SList, ns = ns, TList = TList,
                             lambda = lambda, lambdaFmat = lambdaFmat,
                             verbose = FALSE, ...)
      slh <- c(slh, .LL(Sik, PList[[k]]))
    }
  }
  return(mean(slh))
}



.afcvl <- function(lambda, lambdaFmat, YList, TList, ...) {
  ##############################################################################
  # - (Internal) For at given lambda and lambdaF, compute the approximate
  # - LOOCV loss
  # - lambda > numeric of length 1 giving ridge penalty.
  # - lambdaFmat > numeric matrix giving the fused penalty matrix.
  # - YList   > A list of length K of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - TList   > A list of length K of target matrices the same size
  #             as those of PList. Default is given by default.target.
  # - ...     > Arguments passed to ridgeS.fused
  ##############################################################################

  ns <- sapply(YList, nrow)
  K <- length(ns)
  SList <- lapply(YList, covML)
  PList <- ridgeS.fused(SList = SList, TList = TList, ns = ns,
                       lambda = lambda, lambdaFmat = lambdaFmat,
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



.parseLambda <- function(lambdaFmat) {
  ##############################################################################
  # - A function to parse a character matrix that defines the class of penalty
  #   graphs and unique parameters for cross validation. Returns a list of
  #   indices for each level to be penalized equally.
  #   This list is to be used to construct numeric matrices of penalties.
  # - lambdaFmat > A square K by K character matrix defining the class penalty
  #                matrices to use. Entries with NA, "NA", "" (the empty
  #                string), or "0" are used specify that the pair should be
  #                omitted.
  ##############################################################################

  stopifnot(is.character(lambdaFmat))
  stopifnot(is.matrix(lambdaFmat))
  stopifnot(nrow(lambdaFmat) == ncol(lambdaFmat))

  lambdaFmat[is.na(lambdaFmat)] <- ""
  lambdaFmat[lambdaFmat %in% c("0", "NA")] <- ""

  lvls <- unique(c(lambdaFmat))
  lvls <- lvls[lvls != ""]

  # For each non-empty level get boolean matrices
  parsedLambda <-
    lapply(lvls, function(lvl) which(lvl == lambdaFmat, arr.ind = TRUE))
  names(parsedLambda) <- lvls

  return(parsedLambda)
}



.reconstructLambda <- function(lambdas, parsedLambda, K) {
  ##############################################################################
  # - Reconstruct the numeric penalty matrix lambdaFmat from vector (lambdas)
  #   using output from .parseLambda output.
  # - lambdas      > A numeric vector where the first entry is the
  #                  ridge penalty and the remaining are the
  #                  variable fused penalties.
  # - parsedLambda > A list of length K of matrix indicies.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  get.num <- suppressWarnings(as.numeric(names(parsedLambda)))

  if (length(lambdas) != length(parsedLambda[is.na(get.num)]) + 1) {
    stop("The number of lambdas does not correspond with the number of",
         " non-fixed penalties given i parsedLambda")
  }

  lambdaFmat <- matrix(0, K, K)
  j <- 1
  for (i in seq_along(parsedLambda)) {
    if (is.na(get.num[i])) {
      lambdaFmat[parsedLambda[[i]]] <- lambdas[-1][j]
      j <- j + 1
    } else {
      lambdaFmat[parsedLambda[[i]]] <- get.num[i]
    }
  }
  return(lambdaFmat)
}



optPenalty.fused.LOOCV <- function(YList,
                                   TList,
                                   lambdaMin, lambdaMax,
                                   step1 = 20,
                                   lambdaFMin = lambdaMin,
                                   lambdaFMax = lambdaMax,
                                   step2 = step1,
                                   approximate = FALSE,
                                   verbose = TRUE,
                                   ...) {
  ##############################################################################
  #   Simple (approximate) leave one-out cross validation for the fused ridge
  #   estimator on a grid to determine optimal lambda and lambdaF.
  #   The complete penalty graph is assumed.
  # - YList       > A list of length K of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - TList       > A list of length K of target matrices the same size
  #                 as those of PList. Default is given by default.target.
  # - lambdaMin  > Start of lambda value, the ridge penalty
  # - lambdaMax  > End of lambda value
  # - step1       > Number of evaluations
  # - lambdaFMin  > As lambdaMin for the fused penalty. Default is lambdaMin.
  # - lambdaFMax  > As lambdaMax for the fused penalty. Default is lambdaMax.
  # - step2       > As step1 for the fused penalty. Default is step1.
  # - approximate > Should approximate LOOCV be used? Defaults is FALSE.
  #                 Approximate LOOCV is much faster.
  # - ...         > Arguments passed to ridgeS.fused
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
  lambdas  <- exp(seq(log(lambdaMin),  log(lambdaMax),  length.out = step1))
  lambdaFs <- exp(seq(log(lambdaFMin), log(lambdaFMax), length.out = step2))
  stopifnot(all(is.finite(lambdas)))
  stopifnot(all(is.finite(lambdaFs)))

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
  for (l1 in seq_along(lambdas)) {
    for (l2 in seq_along(lambdaFs)) {
      slh[l1, l2] <-
        cvl(lambda = lambdas[l1], lambdaFmat = matrix(lambdaFs[l2], K, K),
            YList = YList, TList = TList, ...)
      if (verbose){
        cat(sprintf("lambda = %.3f (%d), lambdaF = %.3f (%d), -ll = %.3f\n",
                   lambdas[l1],  l1, lambdaFs[l2], l2, slh[l1, l2]))
      }
    }
  }
  return(list(lambda = lambdas, lambdaF = lambdaFs, fcvl = slh))
}



optPenalty.fused.LOOCVauto <- function(YList,
                                       TList,
                                       lambdaFmat,
                                       approximate = FALSE,
                                       verbose = TRUE,
                                       maxit.ridgeS.fused = 1000,
                                       optimizer = "optim",
                                       maxit.optimizer = 1000,
                                       debug = FALSE,
                                       ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   leave-one-out cross-validation using multi-dimensional optimization
  #   routines.
  #
  # - YList       > A list of length K of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - TList       > A list of length K of target matrices the same size
  #                 as those of PList. Default is given by default.target.
  # - lambdaFmat  > A K by K character matrix defining the class of penalty
  #                 graph to use. The unique elements of lambdaFmat specify the
  #                 penalties to determine. Pairs can be left out using either
  #                 of "", NA, "NA" or "0".
  # - approximate > logical. Should approximate LOOCV be used?
  # - verbose     > logical. Should the function print extra info. Defaults to
  #                 TRUE.
  # - maxit.ridgeS.fused > integer. Max. number of iterations for ridgeS.fused
  # - optimizer          > character giving the stadard optimizer.
  #                        Either "optim" or "nlm".
  # - maxit.optimizer    > integer. Max. number of iterations for the optimizer.
  # - debug              > logical. If TRUE the raw output from the optimizer is
  #                        appended as an attribute to the output.
  # - ...                > arguments passed to the optimizer.
  #
  # The function returns a list of length 4 with entries (1) lambda,
  # (2) lambdaF, (3) the optimal penalty matrix lambdaFmat, and (4) the value of
  # the loss in the optimum. If lambdaFmat is the complete graph, then lambdaF
  # is given. Otherwise lambdaF is NA.
  ##############################################################################

  K <- length(YList)

  # Handle lambdaFmat
  if (missing(lambdaFmat)) {
    lambdaFmat <- matrix("A", K, K)
    diag(lambdaFmat) <- ""
  }

  parsedLambda <- .parseLambda(lambdaFmat)
  n.lambdas <- sum(is.na(suppressWarnings(as.numeric(names(parsedLambda))))) + 1

  # Determine what loss function to use
  # lambdas[1] is the regular ridge penalty, while the remaning lambdas[-1]
  # correspond to the penalty matrix.
  # We also reparameterize to work on log-scale
  if (approximate) {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .afcvl(lambda = elambdas[1],
             lambdaFmat = .reconstructLambda(elambdas, parsedLambda, K),
             YList = YList, TList = TList, maxit = maxit.ridgeS.fused, ...)
    }
  } else {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .fcvl(lambda = elambdas[1],
            lambdaFmat = .reconstructLambda(elambdas, parsedLambda, K),
            YList = YList, TList = TList, maxit = maxit.ridgeS.fused, ...)
    }
  }

  # Get sensible starting value for lambda (choosing lambdaFmat to be zero)
  st <- optimize(function(x) cvl(c(x, rep(0, n.lambdas - 1))),
                 lower = -30, upper = 30)

  # Start lambdaFmat at 0
  lambdas.init <- c(st$minimum, rep(0, n.lambdas - 1))

  if (optimizer == "optim") {

    ans <- optim(lambdas.init, fn = cvl, ...,
                 control = list(trace = verbose, maxit = maxit.optimizer))
    par <- ans$par
    val <- ans$value

  } else if (optimizer == "nlm") {

    ans <- nlm(cvl, lambdas.init, iterlim = maxit.optimizer, ...)
    par <- ans$estimate
    val <- ans$minimum

  }

  # Format optimal values
  opt.lambdas <- exp(par)
  opt.lambdaFmat <- .reconstructLambda(opt.lambdas, parsedLambda, K)
  res <- list(lambda = opt.lambdas[1],
              lambdaF = NA,
              lambdaFmat = opt.lambdaFmat,
              value = val)

  tmp <- unique(opt.lambdaFmat[lower.tri(opt.lambdaFmat)])
  if (length(tmp) == 1) {
    res$lambdaF <- tmp
  }

  if (debug) {
    attr(res, "optim.debug") <- ans
  }

  return(res)
}


##
##
## Help in choosing targets
##
##

default.target.fused <- function(SList, ns, type = "DAIE", equal = TRUE, ...) {
  ##############################################################################
  # Generate a list of (data-driven) targets to use in fused ridge estimation
  # A nice wrapper for default.target
  # - SList > A list of covariance matrices
  # - ns    > A numeric vector of sample sizes corresponding to SList
  # - type  > A character giving the choice of target. See default.target.
  # - equal > logical. If TRUE, all entries in the list are identical and
  #           computed from the pooled estimate. If FALSE, the target is
  #           calculated from each entry in SList.
  # - ...   > Arguments passed to default.target.
  # See also default.target
  ##############################################################################

  stopifnot(is.list(SList))
  stopifnot(is.numeric(ns))
  stopifnot(length(SList) == length(ns))

  if (equal) {
    pooled <- pooledS(SList, ns)
    Tpool <- default.target(pooled, type = type, ...)
    TList <- replicate(length(SList), Tpool, simplify = FALSE)
  } else {
    TList <- lapply(SList, default.target, type = type, ...)
  }
  return(TList)
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
  # Construct the Cartesian product graph from two "character" matrices.
  # - A     > A character matrix where "" signify non-adjacency.
  # - B     > A character matrix where "" signify non-adjacency.
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "X")
  # B <- rags2ridges:::.charAdjMat(factor(letters[4:5]), name = "Y")
  # rags2ridges:::.cartesianProd(A, B)
  ##############################################################################

  AI <- kronecker(.char2num(A), diag(nrow(B)), make.dimnames = TRUE)
  IB <- kronecker(diag(nrow(A)), .char2num(B), make.dimnames = TRUE)
  gprod <- AI + IB  # Kronecker sum

  ans <- kronecker(A, B, FUN = paste0, make.dimnames = TRUE)
  ans[!as.logical(gprod)] <- ""
  return(ans)
}



.tensorProd <- function(A, B) {
  ##############################################################################
  # Construct the Tensor (or categorical) product graph from two "character"
  # matrices.
  # - A     > A character matrix where "" signify non-adjacency.
  # - B     > A character matrix where "" signify non-adjacency.
  # Examples:
  # A <- rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "X")
  # B <- rags2ridges:::.charAdjMat(factor(letters[4:5]), name = "Y")
  # rags2ridges:::.tensorProd(A, B)
  ##############################################################################

  gprod <- kronecker(.char2num(A), .char2num(B), make.dimnames = TRUE)

  ans <- kronecker(A, B, FUN = paste0, make.dimnames = TRUE)
  ans[!as.logical(gprod)] <- ""
  return(ans)
}



default.penalty <- function(K, df,
                            type = c("Complete", "CartesianEqual",
                                     "CartesianUnequal", "TensorProd")) {
  ##############################################################################
  # Select a one of standard penalty matrix types from a dataframe
  # - K     > The number of classes. Can also be list of length K such as
  #           the usual argument "SList".
  #           Can be omitted if 'df' is given.
  # - df    > A data.frame with K rows and with factors in the columns.
  #           Columns of type character are coerced to factors.
  #           Can be omitted when 'type == "Complete"'.
  # - type  > A character giving the type of fused penalty graph to construct.
  #           Should be one of 'Complete' (default), 'CartesianEqual',
  #           'CartesianUnequal', or 'TensorProd' or an unique abbreviation
  #           hereof.
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

  } else if (type == "TensorProd") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.tensorProd, adj.mats)
    return(M)

  } else {

    stop("type =", type, "not implemented yet!")

  }
}






################################################################################
## -------------------------------------------------------------------------- ##
##
## Sparsification and network stats
##
## -------------------------------------------------------------------------- ##
################################################################################


sparsify.fused <- function(PList, ...) {
  ##############################################################################
  # Simple wrapper for sparsify. See help(sparsify).
  # - PList > A list of precision matrices.
  # - ...   > Arguments passed to sparsify.
  ##############################################################################

  return(lapply(PList, sparsify, ...))
}



GGMnetworkStats.fused <- function(PList) {
  ##############################################################################
  # Simple wrapper for GGMnetworkStats. See help(GGMnetworkStats).
  # - PList > A list of sparse precision matrices.
  ##############################################################################

  res <- lapply(PList, GGMnetworkStats, as.table = TRUE)
  if (is.null(names(res))) {
    names(res) <- seq_along(PList)
  }
  return(as.data.frame(res))
}



GGMpathStats.fused <- function(sparsePList, ...) {
  ##############################################################################
  # A wrapper for GGMpathStats in the fused case. See GMMpathStats.
  # - sparsePList > A list of sparsified precision matrices
  # - ...         > Arguments passed to GGMpathStats
  ##############################################################################

  # See if verbose is in ... and set to GGMpathStats default if not
  args <- list(...)
  if (is.null(args[["verbose"]])) {
    verbose <- formals(GGMpathStats)$verbose
  }

  # Run through each class
  res <- vector("list", length(sparsePList))
  names(res) <- names(sparsePList)
  for (k in seq_along(res)) {
    if (verbose) {
      cat("\n\n========================================\n",
          "Class: ", names(res)[k], "\n", sep = "")
    }
    res[[k]] <- GGMpathStats(sparsePList[[k]], ...)
  }
  return(res)
}




