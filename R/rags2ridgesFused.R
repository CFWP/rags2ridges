################################################################################
################################################################################
################################################################################
##
## Fused Ridge estimation with supporting functions for high-dimensional
## precision matrices
##
## Authors: Anders E. Bilgrau, Carel F.W. Peeters, Wessel N. van Wieringen
## Email:	  cf.peeters@vumc.nl
##
################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Section: Support Functions
##
##------------------------------------------------------------------------------
################################################################################
################################################################################


isSymmetricPD <- function(M) {
  ##############################################################################
  # Test if matrix is symmetric postive definite
  # - M > A numeric matrix
  # Returns TRUE if the matrix is symmetric positive definite and FALSE if not.
  ##############################################################################

  nm <- deparse(substitute(M))
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(nm, " is not a numeric matrix")
  }
  if (!isSymmetric(M)) {
    stop(nm, " is not a symmetric matrix")
  }

  chM <- try(chol(M), silent = TRUE)  # M is P.D. iff it has a Cholesky decomp.
  if (inherits(chM, "try-error")) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}

is.Xlist <- function(Xlist, Ylist = FALSE) {
  ##############################################################################
  # Test if generic fused list arguments (such as Slist, Tlist, Plist)
  # are properly formatted
  # - Xlist > A list of covariance matrices or matrices.
  # Returns TRUE if all tests are passed, throws error if not.
  ##############################################################################

  xlist <- deparse(substitute(Xlist))
  if (!is.list(Xlist)) {
    stop(xlist, "should be a list")
  }
  if (!all(sapply(Xlist, is.matrix))) {
    stop("All elements of ", xlist, " should be matrices")
  }
  if (!all(sapply(Xlist, is.numeric))) {
    stop("All elements of ", xlist, " should be numeric matrices")
  }
  if (length(unique(c(sapply(Xlist, dim)))) != 1L) {
    stop("All matrices in ", xlist,
         " should be square and have the same size.")
  }
  if (!all(sapply(Xlist, isSymmetricPD))) {
    stop("All matrices in ", xlist, " should be symmetric and positive ",
         "definite.")
  }
  if (!all(sapply(seq_along(Xlist),
                  function(i) identical(dimnames(Xlist[[1]]),
                                        dimnames(Xlist[[i]]))))) {
    stop("dimnames of the elements of ", xlist, " are not identical")
  }

  return(TRUE)
}


# is.Xlist(dlbcl.S)
# is.Xlist(dlbcl.T)
# length(dlbcl.ns) == length(dlbcl.S)
# length(dlbcl.ns) == length(dlbcl.T)
# all(dim(optimal.penalties$lambdaFmat) == length(dlbcl.ns))



default.target.fused <- function(Slist, ns, type = "DAIE", equal = TRUE, ...) {
  ##############################################################################
  # Generate a list of (data-driven) targets to use in fused ridge estimation
  # A nice wrapper for default.target
  # - Slist > A list of covariance matrices
  # - ns    > A numeric vector of sample sizes corresponding to Slist
  # - type  > A character giving the choice of target. See default.target.
  # - equal > logical. If TRUE, all entries in the list are identical and
  #           computed from the pooled estimate. If FALSE, the target is
  #           calculated from each entry in Slist.
  # - ...   > Arguments passed to default.target.
  # See also default.target
  ##############################################################################

  stopifnot(is.list(Slist))
  stopifnot(is.numeric(ns))
  stopifnot(length(Slist) == length(ns))

  if (equal) {
    pooled <- pooledS(Slist, ns)
    Tpool <- default.target(pooled, type = type, ...)
    Tlist <- replicate(length(Slist), Tpool, simplify = FALSE)
  } else {
    Tlist <- lapply(Slist, default.target, type = type, ...)
  }
  return(Tlist)
}



createS <- function(n, p,
                    topology = "identity",
                    dataset = FALSE,
                    precision = FALSE,
                    nonzero = 0.25,
                    m = 1L,
                    banded.n = 2L,
                    invwishart = FALSE,
                    nu = p + 1,
                    Plist) {
  ##############################################################################
  # - Simulate one or more random symmetric square (or datasets) matrices from
  #   various models.
  # - n          > A vector of sample sizes
  # - p          > An integer giving the dimension. p should be greater than
  #                or equal to 2.
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
  # - invwishart > logical. If TRUE the constructed precision matrix is
  #                used as the scale matrix of an inverse Wishart distribution
  #                and class covariance matrices are drawn from this
  #                distribution.
  # - nu         > The degrees of freedom in the inverse wishart distribution.
  #                A large nu implies high class homogeneity.
  #                Must be greater than p + 1.
  # - Plist      > A list of user-supplied precision matrices. Should be the
  #                same length as n.
  # - Returns a list of matrices if length(n) > 1. The output is simplified if
  #   n has length 1 where only the matrix is returned
  ##############################################################################

  if (missing(p) && !missing(Plist)) {
    p <- nrow(Plist[[1]])
  }
  stopifnot(p > 1)
  stopifnot(m >= 1)
  G <- length(n)

  if (dataset && precision) {
    stop("dataset and precision cannot be TRUE at the same time.")
  }

  if (invwishart && missing(nu)) {
    stop("argument 'nu' is missing. Supply the degrees of freedom 'nu' for ",
         "the inverse Wishart distribution.")
  }

  topology <- match.arg(topology,
                        c("identity", "star", "clique", "complete",
                          "chain", "banded", "Barabassi", "small-world",
                          "scale-free", "Watts-Strogatz", "random-graph",
                          "Erdos-Renyi"))

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
  ans <- vector("list", G)
  names(ans) <- paste0("class", seq_len(G))
  for (g in seq_len(G)) {

    if (!missing(Plist)) {
      stopifnot(length(Plist) == length(n))
      stopifnot(nrow(Plist[[g]]) == ncol(Plist[[g]]))
      stopifnot(nrow(Plist[[g]]) == p)

      Sg <- solve(Plist[[g]])
    } else if (invwishart) {
      stopifnot(nu - p - 1 > 0)
      Sg <- drop(.armaRInvWishart(n = 1, psi = (nu - p - 1)*S, nu = nu))
    } else {
      Sg <- S
    }

    if (precision) {

      if (invwishart) {
        ans[[g]] <- solve(Sg)
      } else {
        ans[[g]] <- P
      }

    } else {

      ans[[g]] <- rmvnormal(n = n[g], mu = rep(0, p), sigma = Sg)

      if (!dataset) {
        ans[[g]] <- covML(ans[[g]])
      }

    }

    if (p <= 17576) {  # Only give names for "small" dimensions
      colnames(ans[[g]]) <- nms[1:p]
      if (!dataset) {
        rownames(ans[[g]]) <- nms[1:p]
      }
    }
  }

  if (G == 1) {  # Simplify output if ns is length 1
    ans <- ans[[1]]
  }

  return(ans)
}



pooledS <- function(Slist, ns, mle = TRUE) {
  ##############################################################################
  # - Computes the pooled covariance estimate
  # - Slist > A list sample covariance matrices for each class
  # - ns    > A vector of sample sizes of the same length as Slist.
  # - mle   > logical. If TRUE the biased MLE is used. If FALSE, the biased
  #           corrected estimate is used.
  ##############################################################################
  ans <- .armaPooledS(Slist = Slist, ns = ns, mle = as.numeric(mle))
  dimnames(ans) <- dimnames(Slist[[1]])
  return(ans)
}



KLdiv.fused <- function(MtestList, MrefList, StestList, SrefList, ns,
                        symmetric = FALSE) {
  ##############################################################################
  # - Function that calculates the "weigthed" Kullback-Leibler divergence
  #   between multiple paired normal distributions
  # - MtestList > A list of mean vectors approximating. Assumed to be zero
  #               vectors if not supplied.
  # - MrefList  > A list of mean vectors 'true'/reference. Assumed to be zero
  #               vectors if not supplied.
  # - StestList > A list of covariance matrix approximating
  # - SrefList  > A list of covariance matrix 'true'/reference
  # - symmetric > logical indicating if original symmetric version of KL div.
  #               should be calculated.
  ##############################################################################

  if (missing(MtestList)) {
    MtestList <- replicate(length(StestList), rep(0, nrow(StestList[[1]])),
                           simplify = FALSE)
  }
  if (missing(MrefList)) {
    MrefList <- replicate(length(StestList), rep(0, nrow(StestList[[1]])),
                          simplify = FALSE)
  }
  KLdivs <- mapply(KLdiv, MtestList, MrefList, StestList, SrefList,
                   MoreArgs = list(symmetric = symmetric))

  return(sum(ns*KLdivs)/sum(ns))
}



.FLL <- function(Slist, Plist, ns){
  ##############################################################################
  # - Function that computes the value of the (negative) combined log-likelihood
  # - Slist > A list sample covariance matrices for each class
  # - Plist > A list of the same length as (Slist) of precision matrices
  #   (possibly regularized inverse of covariance or correlation matrices)
  # - ns > A vector of sample sizes of the same length as Slist.
  ##############################################################################

  LLs <- mapply(.LL, Slist, Plist)
  return(sum(ns*LLs))
}



.PFLL <- function(Slist, Plist, ns, Tlist, lambda, lambdaFmat){
  ##############################################################################
  # - Function that computes the value of the (negative) penalized combined
  #   log-likelihood
  # - Slist   > A list sample covariance matrices for each class
  # - Plist   > A list of the same length as (Slist) of precision matrices
  #             (possibly regularized inverse of covariance  matrices)
  # - ns      > A vector of sample sizes of the same length as Slist.
  # - Tlist   > list of target matrices
  # - lambda  > ridge penalty
  # - lambdaFmat > fused penalty matrix
  ##############################################################################

  penalty <- 0
  for (g1 in seq_along(Slist)) {
    for (g2 in seq_len(g1)) {
      if (g1 == g2) { # Ridge penalty
        penalty <- penalty +
          lambda*.FrobeniusLoss(Slist[[g1]], Tlist[[g1]])
      } else {  # Fused contribution
        penalty <- penalty +
          lambdaFmat[g1, g2]*.FrobeniusLoss(Slist[[g1]] - Tlist[[g1]],
                                            Slist[[g2]] - Tlist[[g2]])
      }
    }
  }
  penalty <- penalty/2

  ans <- .FLL(Slist, Plist, ns) + penalty
  return(ans)
}



################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Section: The fused ridge estimator
##
##------------------------------------------------------------------------------
################################################################################
################################################################################



.fusedUpdateI <- function(g0, Plist, Slist, Tlist, ns, lambda, lambdaFmat) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate. The scheme I approach.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda > The ridge penalty (a postive number).
  # - lambdaFmat > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambdaFmat[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: The update function seems to work ok for large lambdaFmat.
  #   However, for very large lambdaFmat (> 1e154) the exception that the
  #   .armaRidgeP returns the target because of an exception. Which is wrong
  #   in the fused case.
  ##############################################################################

  diag(lambdaFmat) <- 0  # Make sure there's zeros in the diagonal
  a <- (sum(lambdaFmat[g0, ]) + lambda)/ns[g0]
  b <- lambdaFmat[g0, -g0]/ns[g0]

  OmT <- mapply(`-`, Plist[-g0], Tlist[-g0], SIMPLIFY = FALSE) # Omega - Target
  OmT <- mapply(`*`, b, OmT, SIMPLIFY = FALSE)
  S0 <- Slist[[g0]] - Reduce(`+`, OmT)
  return(.armaRidgeP(S0, target = Tlist[[g0]], lambda = a))
}



.fusedUpdateII <- function(g0, Plist, Slist, Tlist, ns, lambda, lambdaFmat) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate -- using the alternative II update scheme.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda > The ridge penalty (a postive number).
  # - lambdaFmat > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambdaFmat[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: This update seems to work very poorly for large lambdaFmat
  ##############################################################################

  p <- nrow(Plist[[1]])
  lambdaa <- (lambda + sum(lambdaFmat[g0, -g0]))/ns[g0]
  b <- (lambda + sum(lambdaFmat[g0, -g0]) - 1)/ns[g0]

  Psum <- Tsum <- matrix(0, p, p)
  for (g in setdiff(seq_along(Plist), g0)) {
    Psum <- Psum + lambdaFmat[g0, g]*Plist[[g]]
    Tsum <- Tsum + (lambdaFmat[g0, g]/ns[g0])*Tlist[[g]]
  }
  Sbar <- Slist[[g0]] + b*Psum + Tsum
  Tbar <- Tlist[[g0]] + Psum
  return(.armaRidgeP(Sbar, target = Tbar, lambda = lambdaa))
}



.fusedUpdateIII <- function(g0, Plist, Slist, Tlist, ns, lambda, lambdaFmat) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate -- using the alternative III update scheme.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda > The ridge penalty (a postive number).
  # - lambdaFmat > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambdaFmat[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: This update function seems to work very well for large lambdaFmat.
  #   For very large lambdaFmat (> 1e154) the exception triggered in the
  #   .armaRidgeP returns the target because of an exception. However, in this
  #   updating scheme, that is also correct.
  ##############################################################################

  lambdasum <- (lambda + sum(lambdaFmat[g0, -g0]))
  lambdaa <- lambdasum/ns[g0]

  Tbar <- Tlist[[g0]]
  for (g in setdiff(seq_along(Plist), g0)) {
    Tbar <- Tbar + (lambdaFmat[g0, g]/lambdasum)*(Plist[[g]] - Tlist[[g]])
  }

  return(.armaRidgeP(Slist[[g0]], target = Tbar, lambda = lambdaa))
}



ridgeP.fused <- function(Slist, ns, Tlist = default.target.fused(Slist, ns),
                         lambda, lambdaFmat, lambdaF, Plist,
                         maxit = 100L, verbose = TRUE, eps = 1e-4) {
  ##############################################################################
  # - The user function for the fused ridge estimate for a given
  #   lambda and lambdaFmat.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda > The ridge penalty (a postive number)
  # - lambdaFmat > The G by G symmetric adjacency matrix fused penalty graph
  #             with non-negative entries where lambdaFmat[g1, g2] determine the
  #             retainment of similarities between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  # - lambdaF > The non-negative fused penalty. Alternative to lambdaFmat if
  #             all pairwise penalties equal.
  # - Plist   > A list of length G giving the initial estimates. If not supplied
  #             the ridge of the pooled estimate is used.
  # - maxit   > integer. The maximum number of interations, default is 100.
  # - verbose > logical. Should the function print extra info. Defaults to TRUE.
  # - eps     > numeric. A positive convergence criterion. Default is 1e-4.
  ##############################################################################

  stopifnot(length(Slist) == length(Tlist))
  G <- length(Slist)  # Number of groups

  # Initialize estimates with the regular ridges from the pooled covariance
  if (missing(Plist)) {
    Spool <- pooledS(Slist, ns, mle = FALSE)
    Plist <- list()
    for (i in seq_len(G)) {
      Plist[[i]] <- .armaRidgeP(Spool, target = Tlist[[i]],
                               lambda = G*lambda/sum(ns))
    }
  }
  stopifnot(length(Slist) == length(Plist))

  if (!missing(lambdaFmat) && !missing(lambdaF)) {
    stop("Supply only either lambdaFmat or lambdaF.")
  } else if (missing(lambdaFmat) && missing(lambdaF)) {
    stop("Either lambdaFmat or lambdaF must be given.")
  } else if (missing(lambdaFmat) && !missing(lambdaF)) {
    lambdaFmat <- matrix(lambdaF, G, G)
  }

  # Overwrite the starting estimate with the fused estimate
  Plist <-
    .armaRidgeP_fused(Slist = Slist, ns = ns, Tlist = Tlist, lambda = lambda,
                      lambdaFmat = lambdaFmat,Plist = Plist, maxit = maxit,
                      eps = eps, verbose = verbose)

  if (i == maxit + 1) {
    warning("Maximum iterations (", maxit, ") hit")
  }

  # Keep dimnames and names
  for (g in seq_along(Slist)) {
    dimnames(Plist[[g]]) <- dimnames(Slist[[g]])
  }
  names(Plist) <- names(Slist)


  return(Plist)
}



################################################################################
################################################################################
## ----------------------------------------------------------------------------
##
## Section: LOOCV (and approximation) in the fused setting
##
## -----------------------------------------------------------------------------
################################################################################
################################################################################


.fcvl <- function(lambda, lambdaFmat, Ylist, Tlist, ...) {
  ##############################################################################
  # - A function to perform fused LOOCV
  # - lambda > numeric of length 1 giving ridge penalty.
  # - lambdaFmat > numeric matrix giving the fused penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ...     > Arguments passed to ridgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  Slist.org <- lapply(Ylist, covML)

  for (g in seq_len(G)) {
    ng <- nrow(Ylist[[g]])
    slh <- numeric()
    for (i in seq_len(ng)) {
      ns <- ns.org
      ns[g] <- ns[g] - 1
      Slist <- Slist.org
      Slist[[g]] <- covML(Ylist[[g]][-i, , drop = FALSE])
      Sig    <- crossprod(Ylist[[g]][i,  , drop = FALSE])
      Plist  <- ridgeP.fused(Slist = Slist, ns = ns, Tlist = Tlist,
                             lambda = lambda, lambdaFmat = lambdaFmat,
                             verbose = FALSE, ...)
      slh <- c(slh, .LL(Sig, Plist[[g]]))
    }
  }
  return(mean(slh))
}



.afcvl <- function(lambda, lambdaFmat, Ylist, Tlist, ...) {
  ##############################################################################
  # - (Internal) For at given lambda and lambdaF, compute the approximate
  # - LOOCV loss
  # - lambda > numeric of length 1 giving ridge penalty.
  # - lambdaFmat > numeric matrix giving the fused penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ...     > Arguments passed to ridgeP.fused
  ##############################################################################

  ns <- sapply(Ylist, nrow)
  G <- length(ns)
  Slist <- lapply(Ylist, covML)
  Plist <- ridgeP.fused(Slist = Slist, Tlist = Tlist, ns = ns,
                       lambda = lambda, lambdaFmat = lambdaFmat,
                       verbose = FALSE, ...)
  n.tot <- sum(ns)
  nll <- .FLL(Slist = Slist, Plist = Plist, ns)/n.tot
  denom <- n.tot*(n.tot - 1)
  bias <- 0
  for (g in seq_along(ns)) {
    for (i in seq_len(ns[g])) {
      Sig <- crossprod(Ylist[[g]][i, , drop = FALSE])
      fac1 <- diag(nrow(Sig)) - Sig %*% Plist[[g]]
      fac2 <- Plist[[g]] %*% (Slist[[g]] - Sig) %*% Plist[[g]]
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
  # - lambdaFmat > A square G by G character matrix defining the class penalty
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



.reconstructLambda <- function(lambdas, parsedLambda, G) {
  ##############################################################################
  # - Reconstruct the numeric penalty matrix lambdaFmat from vector (lambdas)
  #   using output from .parseLambda output.
  # - lambdas      > A numeric vector where the first entry is the
  #                  ridge penalty and the remaining are the
  #                  variable fused penalties.
  # - parsedLambda > A list of length G of matrix indicies.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  get.num <- suppressWarnings(as.numeric(names(parsedLambda)))

  if (length(lambdas) != length(parsedLambda[is.na(get.num)]) + 1) {
    stop("The number of lambdas does not correspond with the number of",
         " non-fixed penalties given i parsedLambda")
  }

  lambdaFmat <- matrix(0, G, G)
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



optPenalty.fused.LOOCVgrid <- function(Ylist,
                                       Tlist,
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
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambdaMin  > Start of lambda value, the ridge penalty
  # - lambdaMax  > End of lambda value
  # - step1       > Number of evaluations
  # - lambdaFMin  > As lambdaMin for the fused penalty. Default is lambdaMin.
  # - lambdaFMax  > As lambdaMax for the fused penalty. Default is lambdaMax.
  # - step2       > As step1 for the fused penalty. Default is step1.
  # - approximate > Should approximate LOOCV be used? Defaults is FALSE.
  #                 Approximate LOOCV is much faster.
  # - ...         > Arguments passed to ridgeP.fused
  # - verbose     > logical. Print extra information. Defaults is TRUE.
  #
  # The function evaluates the loss on a log-equidistant grid.
  ##############################################################################

  if (missing(Tlist)) {  # If Tlist is not provided
    Tlist <- lapply(Ylist, function(Y) default.target(covML(Y)))
  }

  G <- length(Ylist)
  ns <- sapply(Ylist, nrow)

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
  total.n <- sum(sapply(Ylist, nrow))
  for (l1 in seq_along(lambdas)) {
    for (l2 in seq_along(lambdaFs)) {
      slh[l1, l2] <-
        cvl(lambda = lambdas[l1], lambdaFmat = matrix(lambdaFs[l2], G, G),
            Ylist = Ylist, Tlist = Tlist, ...)
      if (verbose){
        cat(sprintf("lambda = %.3f (%d), lambdaF = %.3f (%d), -ll = %.3f\n",
                   lambdas[l1],  l1, lambdaFs[l2], l2, slh[l1, l2]))
      }
    }
  }
  return(list(lambda = lambdas, lambdaF = lambdaFs, fcvl = slh))
}



optPenalty.fused.LOOCVauto <- function(Ylist,
                                       Tlist,
                                       lambdaFmat,
                                       approximate = FALSE,
                                       verbose = TRUE,
                                       maxit.ridgeP.fused = 1000,
                                       optimizer = "optim",
                                       maxit.optimizer = 1000,
                                       debug = FALSE,
                                       ...) {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   leave-one-out cross-validation using multi-dimensional optimization
  #   routines.
  #
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambdaFmat  > A G by G character matrix defining the class of penalty
  #                 graph to use. The unique elements of lambdaFmat specify the
  #                 penalties to determine. Pairs can be left out using either
  #                 of "", NA, "NA" or "0".
  # - approximate > logical. Should approximate LOOCV be used?
  # - verbose     > logical. Should the function print extra info. Defaults to
  #                 TRUE.
  # - maxit.ridgeP.fused > integer. Max. number of iterations for ridgeP.fused
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

  G <- length(Ylist)

  # Handle lambdaFmat
  if (missing(lambdaFmat)) {
    lambdaFmat <- matrix("A", G, G)
    diag(lambdaFmat) <- ""
  }

  parsedLambda <- .parseLambda(lambdaFmat)
  suppressWarnings({
    n.lambdas <- sum(is.na(as.numeric(names(parsedLambda)))) + 1
  })

  # Determine what loss function to use
  # lambdas[1] is the regular ridge penalty, while the remaning lambdas[-1]
  # correspond to the penalty matrix.
  # We also reparameterize to work on log-scale
  if (approximate) {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .afcvl(lambda = elambdas[1],
             lambdaFmat = .reconstructLambda(elambdas, parsedLambda, G),
             Ylist = Ylist, Tlist = Tlist, maxit = maxit.ridgeP.fused, ...)
    }
  } else {
    cvl <- function(lambdas, ...) {
      elambdas <- exp(lambdas)
      .fcvl(lambda = elambdas[1],
            lambdaFmat = .reconstructLambda(elambdas, parsedLambda, G),
            Ylist = Ylist, Tlist = Tlist, maxit = maxit.ridgeP.fused, ...)
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
  opt.lambdaFmat <- .reconstructLambda(opt.lambdas, parsedLambda, G)
  dimnames(opt.lambdaFmat) <- dimnames(lambdaFmat)

  # Compute estimate at optimal values
  optPlist <- ridgeP.fused(Slist = lapply(Ylist, covML), ns = ns,
                           Tlist = Tlist, lambda = opt.lambdas,
                           lambdaFmat = opt.lambdaFmat,
                           maxit = maxit.ridgeP.fused, verbose = FALSE)

  # Construct results
  res <- list(lambda = opt.lambdas[1],
              lambdaF = NA,
              lambdaFmat = opt.lambdaFmat,
              Plist = optPlist,
              value = val)

  res$lambdaF <- unique(c(opt.lambdaFmat))
  names(res$lambdaF) <- lambdaFmat[match(res$lambdaF, opt.lambdaFmat)]

  if (debug) {
    attr(res, "optim.debug") <- ans
  }

  return(res)
}



optPenalty.fused <- function() {
  ##############################################################################
  # - Selection of the optimal penalties w.r.t. to (possibly approximate)
  #   leave-one-out cross-validation using multi-dimensional optimization
  #   routines.
  #
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambdaFmat  > A G by G character matrix defining the class of penalty
  #                 graph to use. The unique elements of lambdaFmat specify the
  #                 penalties to determine. Pairs can be left out using either
  #                 of "", NA, "NA" or "0".
  # - approximate > logical. Should approximate LOOCV be used?
  # - verbose     > logical. Should the function print extra info. Defaults to
  #                 TRUE.
  # - maxit.ridgeP.fused > integer. Max. number of iterations for ridgeP.fused
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

  # A wrapper for
  #  optPenalty.fused.LOOCVauto
  #  optPenalty.fused.LOOCV
  message("to be implemented")
  return(-1)
}



################################################################################
################################################################################
## -----------------------------------------------------------------------------
##
## Section: Automatic penalty matrix constructor
##
## -----------------------------------------------------------------------------
################################################################################
################################################################################


.charAdjMat <- function(fac, name = "X") {
  ##############################################################################
  # Create a complete character adjacency matrix from a factor. This function
  # is used in the constructing the character penalty matrix in default.penalty.
  # - fac  > A factor of some length.
  # - name > A character giving the text which should appear in the adjacent
  #          entries. If not a character, the object name of fac is used.
  # Examples:
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "Y")
  # rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = NULL)
  ##############################################################################

  G <- nlevels(fac)
  if (is.character(name)) {
    lab <- name
  } else {
    lab <- deparse(substitute(fac))
  }
  M <- matrix(lab, G, G)
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



default.penalty <- function(G, df,
                            type = c("Complete", "CartesianEqual",
                                     "CartesianUnequal", "TensorProd")) {
  ##############################################################################
  # Select a one of standard penalty matrix types from a dataframe
  # - G     > The number of classes. Can also be list of length G such as
  #           the usual argument "Slist".
  #           Can be omitted if 'df' is given.
  # - df    > A data.frame with G rows and with factors in the columns.
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

  if (missing(G) && !missing(df)) {
    G <- nrow(df)
  }

  if (is.data.frame(G)) {
    df <- G
    G <- nrow(df)
  }

  if (missing(G) && !missing(df)) {
    G <- nrow(df)
  }

  if (is.list(G)) {
    G <- length(G)
  }

  if (missing(df)) {
    if (type != "Complete") {
      warning("No data.frame 'df' given and 'type' does not equal 'Complete'.",
              " Setting 'type' to 'Complete'")
      type <- "Complete"
    }
    df <- data.frame(Class = factor(seq_len(G)))
  }

  if (!all(sapply(df, is.factor))) {
    stop("Not all columns in the data.frame 'df' are factors")
  }

  stopifnot(G == nrow(df))

  if (type == "Complete") {

    M <- matrix("f", G, G)
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
################################################################################
## -----------------------------------------------------------------------------
##
## Section: Sparsification and network stats
##
## -----------------------------------------------------------------------------
################################################################################
################################################################################


sparsify.fused <- function(Plist, ...) {
  ##############################################################################
  # Simple wrapper for sparsify. See help(sparsify).
  # - Plist > A list of precision matrices.
  # - ...   > Arguments passed to sparsify.
  ##############################################################################

  return(lapply(Plist, sparsify, ...))
}



GGMnetworkStats.fused <- function(Plist) {
  ##############################################################################
  # Simple wrapper for GGMnetworkStats. See help(GGMnetworkStats).
  # - Plist > A list of sparse precision matrices.
  ##############################################################################

  res <- lapply(Plist, GGMnetworkStats, as.table = TRUE)
  if (is.null(names(res))) {
    names(res) <- seq_along(Plist)
  }
  return(as.data.frame(res))
}



GGMpathStats.fused <- function(sparsePlist, ...) {
  ##############################################################################
  # A wrapper for GGMpathStats in the fused case. See GMMpathStats.
  # - sparsePlist > A list of sparsified precision matrices
  # - ...         > Arguments passed to GGMpathStats
  ##############################################################################

  # See if verbose is in ... and set to GGMpathStats default if not
  args <- list(...)
  if (is.null(args[["verbose"]])) {
    verbose <- formals(GGMpathStats)$verbose
  }

  # Run through each class
  res <- vector("list", length(sparsePlist))
  names(res) <- names(sparsePlist)
  for (g in seq_along(res)) {
    if (verbose) {
      cat("\n\n========================================\n",
          "Class: ", names(res)[g], "\n", sep = "")
    }
    res[[g]] <- GGMpathStats(sparsePlist[[g]], ...)
  }
  return(res)
}



################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Section: Miscellaneous
##
##------------------------------------------------------------------------------
################################################################################
################################################################################


.JayZScore <- function() {
  ##############################################################################
  # - The truth
  ##############################################################################

  cat("
      ###################################################
      ###################################################
      I'm from rags to ridges baby I ain't dumb
      I got 99 problems but the ridge ain't one.
      - Jay Z-score
      ###################################################
      ################################################### \n")
}




