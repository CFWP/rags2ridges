################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module B: rags2ridges Fused
##
##------------------------------------------------------------------------------
################################################################################
################################################################################



#' Test for symmetric positive (semi-)definiteness
#'
#' Function to test if a \code{matrix} is symmetric positive (semi)definite or
#' not.
#'
#' Tests positive definiteness by Cholesky decomposition.  Tests positive
#' semi-definiteness by checking if all eigenvalues are larger than
#' \eqn{-\epsilon|\lambda_1|} where \eqn{\epsilon} is the tolerance and
#' \eqn{\lambda_1} is the largest eigenvalue.
#'
#' @param M A square symmetric matrix.
#' @param tol A numeric giving the tolerance for determining positive
#'   semi-definiteness.
#'
#' @return Returns a \code{logical} value. Returns \code{TRUE} if the \code{M}
#'   is symmetric positive (semi)definite and \code{FALSE} if not.  If \code{M}
#'   is not even symmetric, the function throws an error.
#'
#' @author Anders Ellern Bilgrau Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{isSymmetric}}
#'
#' @examples
#' A <- matrix(rnorm(25), 5, 5)
#' \dontrun{
#' isSymmetricPD(A)
#' }
#' B <- symm(A)
#' isSymmetricPD(B)
#'
#' C <- crossprod(B)
#' isSymmetricPD(C)
#'
#' isSymmetricPSD(C)
#'
#' @export
isSymmetricPD <- function(M) {

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



#' @rdname isSymmetricPD
#' @details While \code{isSymmetricPSD} returns \code{TRUE} if the matrix is
#'   symmetric positive definite and \code{FASLE} if not. In practice, it tests
#'   if all eigenvalues are larger than -tol*|l| where l is the largest
#'   eigenvalue. More
#'   \href{https://scicomp.stackexchange.com/questions/12979/testing-if-a-matrix-is-positive-semi-definite}{here.}
#'
#' @export
isSymmetricPSD <- function(M, tol = 1e-4) {
  nm <- deparse(substitute(M))
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(nm, " is not a numeric matrix")
  }
  if (!isSymmetric(M)) {
    stop(nm, " is not a symmetric matrix")
  }

  # TODO: Improve speed
  evals <- eigen(M, symmetric = TRUE)$values # M is PSD iff eigenvals >= 0 SLOW!
  if (all(evals > -tol*abs(max(evals)))) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}



#' Test if fused list-formats are correctly used
#'
#' Function to check if the argument submits to the various \code{list}-formats
#' used by the fused ridge estimator and related functions are correct. That
#' is, it tests if generic fused list arguments (such as \code{Slist},
#' \code{Tlist}, \code{Plist}, \code{Ylist}) are properly formatted.
#'
#' @param Xlist A \code{list} of precision matrices of equal size
#'   (\code{Plist}), sample covariance matrices (\code{Slist}), data matrices
#'   (\code{Ylist})
#' @param Ylist \code{logical}. Is \code{Xlist} a \code{list} of data matrices
#'   with the same number of columns (\code{Ylist}).
#' @param semi \code{logical}. Should the matrices in the list be tested to be
#'   positive semi definite or positive definite?
#'
#' @return Returns \code{TRUE} if all tests are passed, throws error if not.
#'
#' @author Anders Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel N.
#'   van Wieringen
#'
#' @seealso \code{\link{ridgeP.fused}}, \code{\link{optPenalty.fused}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#'   van Wieringen, W.N. & Peeters, C.F.W. (2016).  Ridge Estimation of Inverse
#'   Covariance Matrices from High-Dimensional Data, Computational Statistics &
#'   Data Analysis, vol. 103: 284-303.  Also available as arXiv:1403.0904v3
#'   [stat.ME].
#'
#' @examples
#' Slist <- createS(n = c(4, 6, 9), p = 10)
#' is.Xlist(Slist, semi = TRUE)
#' @export is.Xlist
is.Xlist <- function(Xlist, Ylist = FALSE, semi = FALSE) {

  xlist <- deparse(substitute(Xlist))
  if (!is.list(Xlist)) {
    stop(xlist, " should be a list")
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
  if (semi) {
    if (!all(sapply(Xlist, isSymmetricPSD))) {
      stop("All matrices in ", xlist, " should be symmetric and positive ",
           "semi definite.")
    }
  } else {
    if (!all(sapply(Xlist, isSymmetricPD))) {
      stop("All matrices in ", xlist, " should be symmetric and positive ",
           "definite.")
    }
  }
  if (!all(sapply(seq_along(Xlist),
                  function(i) identical(dimnames(Xlist[[1]]),
                                        dimnames(Xlist[[i]]))))) {
    stop("dimnames of the elements of ", xlist, " are not identical")
  }

  return(TRUE)
}



#' Generate data-driven targets for fused ridge estimation
#'
#' Generates a list of (data-driven) targets to use in fused ridge estimation.
#' Simply a wrapper for \code{\link{default.target}}.
#'
#'
#' @param Slist A \code{list} of length \eqn{K} of \code{numeric} covariance
#'   matrices of the same size for \eqn{K} classes.
#' @param ns A \code{numeric} vector of sample sizes corresponding to the
#'   entries of \code{Slist}.
#' @param type A \code{character} giving the choice of target to construct. See
#'   \code{\link{default.target}} for the available options. Default is
#'   \code{"DAIE"}.
#' @param by A \code{character} vector with the same length as \code{Slist}
#'   specifying which groups should share target.  For each unique entry of
#'   \code{by} a target is constructed.  If omitted, the default is to assign a
#'   unique target to each class.  If not given as a \code{character} coercion
#'   into one is attempted.
#' @param \dots Arguments passed to \code{\link{default.target}}.
#'
#' @return A \code{list} of \eqn{K} covariance target matrices of the same size.
#'
#' @author Anders E. Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel
#'   N. van Wieringen
#'
#' @seealso \code{\link{default.target}}
#'
#' @examples
#' # Make some toy data
#' ns <- c(3, 4)  # Two classes with sample size 3 and 4
#' Slist <- createS(ns, p = 3)  # Generate two 3-dimensional covariance matrices
#' Slist
#'
#' # Different choices:
#' default.target.fused(Slist, ns)
#' default.target.fused(Slist, ns, by = seq_along(Slist)) # The same as before
#' default.target.fused(Slist, ns, type = "Null")
#' default.target.fused(Slist, ns, type = "DAPV")
#' default.target.fused(Slist, ns, type = "DAPV", by = rep(1, length(Slist)))
#'
#'
#' # Make some (more) toy data
#' ns <- c(3, 4, 6, 7)  # Two classes with sample size 3 and 4
#' Slist <- createS(ns, p = 2)  # Generate four 2-dimensional covariance matrices
#'
#' # Use the same target in class 1 and 2, but another in class 3 and 4:
#' default.target.fused(Slist, ns, by = c("A", "A", "B", "B"))
#'
#' @export default.target.fused
default.target.fused <- function(Slist, ns, type = "DAIE", by, ...) {
  stopifnot(is.list(Slist))
  stopifnot(is.numeric(ns))
  stopifnot(length(Slist) == length(ns))

  if (missing(by)) {
    by <- seq_along(Slist)
  }
  stopifnot(length(by) == length(Slist))
  by <- as.character(by)
  stopifnot(!anyNA(by) && is.character(by))

  Tlist <- vector("list", length(Slist))
  names(Tlist) <- names(Slist)
  for (lvl in unique(by)) {
    get <- lvl == by
    pooled <- pooledS(Slist, ns, subset = get)
    Tpool <- default.target(pooled, type = type, ...)
    Tlist[get] <- replicate(sum(get), Tpool, simplify = FALSE)
  }

  return(Tlist)
}









#' Simulate sample covariances or datasets
#'
#' Simulate data from a p-dimensional (zero-mean) gaussian graphical model (GGM)
#' with a specified (or random) topology and return the sample covariance matrix
#' or matrices. Can also return the original simulated data or underlying
#' precision matrix.
#'
#' The data is simulated from a zero-mean \code{p}-dimensional multivariate
#' gaussian distribution with some precision matrix determined by the argument
#' \code{topology} which defines the GGM. If \code{precision} is \code{TRUE} the
#' population precision matrix is returned. This is useful to see what the
#' actual would-be-used precision matrices are. The available values of
#' \code{topology} are described below. Unless otherwise stated the diagonal
#' entries are always one. If \code{m} is 2 or greater block diagonal precision
#' matrices are constructed and used. \itemize{ \item \code{"identity"}: uses
#' the identity matrix (\code{diag(p)}) as precision matrix.  Corresponds to no
#' conditional dependencies.  \item \code{"star"}: simulate from a star
#' topology. Within each block the first node is selected as the "hub". The
#' off-diagonal entries \eqn{(1,j)} and \eqn{(j,1)} values taper off with the
#' value \eqn{1/(j + 1)}.  \item \code{"clique"}: simulate from clique topology
#' where each block is a complete graph with off-diagonal elements equal to
#' \code{nonzero}.  \item \code{"complete"}: alias for (and identical to)
#' \code{"clique"}.  \item \code{"chain"}: simulate from a chain topology where
#' the precision matrix is a tridiagonal matrix with off-diagonal elements (in
#' each block) given by argument \code{nonzero}.  \item \code{"banded"}:
#' precision elements \code{(i,j)} are given by \eqn{1/(|i-j|+1)} if \eqn{|i-j|}
#' is less than or equal to \code{banded.n} and zero otherwise. \item
#' \code{"scale-free"}: The non-zero pattern of each block is generated by a
#' Barabassi random graph. Non-zero off-diagonal values are given by
#' \code{nonzero}.  Gives are very "hubby" network.  \item \code{"Barabassi"}:
#' alias for \code{"scale-free"}.  \item \code{"small-world"}: The non-zero
#' pattern of each block is generated by a 1-dimensional Watts-Strogatz random
#' graph with \code{banded.n} starting neighbors and \eqn{5\%} probability of
#' rewiring.  Non-zero off-diagonal values are given by \code{nonzero}. Gives
#' are very "bandy" network.  \item \code{"Watts-Strogatz"}: alias for
#' \code{"small-world"} \item \code{"random-graph"}: The non-zero pattern of
#' each block is generated by a Erdos-Renyi random graph where each edge is
#' present with probability \eqn{1/p}.  Non-zero off-diagonal values are given
#' by \code{nonzero}.  \item \code{"Erdos-Renyi"}: alias for
#' \code{"random-graph"} } When \code{n} has length greater than 1, the datasets
#' are generated i.i.d. given the topology and number of blocks.
#'
#' Arguments \code{invwishart} and \code{nu} allows for introducing class
#' homogeneity. Large values of \code{nu} imply high class homogeneity.
#' \code{nu} must be greater than \code{p + 1}. More precisely, if
#' \code{invwishart == TRUE} then the constructed precision matrix is used as
#' the scale parameter in an inverse Wishart distribution with \code{nu} degrees
#' of freedom. Each class covariance is distributed according to this inverse
#' Wishart and independent.
#'
#' @param n A \code{numeric} vector giving number of samples. If the length is
#'   larger than 1, the covariance matrices are returned as a list.
#' @param p A \code{numeric} of length 1 giving the dimension of the
#'   samples/covariance.
#' @param topology character. The topology to use for the simulations. See the
#'   details.
#' @param dataset A \code{logical} value specifying whether the sample
#'   covariance or the simulated data itself should be returned.
#' @param precision A \code{logical} value. If \code{TRUE} the constructed
#'   precision matrix is returned.
#' @param nonzero A \code{numeric} of length 1 giving the value of the nonzero
#'   entries used in some topologies.
#' @param m A \code{integer} giving the number of blocks (i.e. conditionally
#'   independent components) to create. If \code{m} is greater than 1, then the
#'   given \code{topology} is used on \code{m} blocks of approximately equal
#'   size.
#' @param banded.n A \code{integer} of length one giving the number of bands.
#'   Only used if \code{topology} is one of \code{"banded"},
#'   \code{"small-world"}, or \code{"Watts-Strogatz"}.
#' @param invwishart \code{logical}. If \code{TRUE} the constructed precision
#'   matrix is used as the scale matrix of an inverse Wishart distribution and
#'   class covariance matrices are drawn from this distribution.
#' @param nu \code{numeric} greater than \code{p + 1} giving the degrees of
#'   freedom in the inverse Wishart distribution.  A large \code{nu} implies
#'   high class homogeneity.  A small \code{nu} near \code{p + 1} implies high
#'   class heterogeneity.
#' @param Plist An optional \code{list} of \code{numeric} matrices giving the
#'   precision matrices to simulate from. Useful when random matrices have
#'   already been generated by setting \code{precision = TRUE}.
#'
#' @return The returned type is dependent on \code{n} and \code{covariance}. The
#'   function generally returns a \code{list} of \code{numeric} matrices with
#'   the same length as \code{n}. If \code{covariance} is \code{FALSE} the
#'   simulated datasets with size \code{n[i]} by \code{p} are given in the
#'   \code{i} entry of the output. If \code{covariance} is \code{TRUE} the
#'   \code{p} by \code{p} sample covariances of the datasets are given. When
#'   \code{n} has length 1 the \code{list} structure is dropped and the matrix
#'   is returned.
#'
#' @author Anders E. Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel
#'   N. van Wieringen
#'
#' @examples
#' ## Generate some simple sample covariance matrices
#' createS(n = 10, p = 3)
#' createS(n = c(3, 4, 5), p = 3)
#' createS(n = c(32, 55), p = 7)
#'
#' ## Generate some datasets and not sample covariance matrices
#' createS(c(3, 4), p = 6, dataset = TRUE)
#'
#' ## Generate sample covariance matrices from other topologies:
#' A <- createS(2000, p = 4, topology = "star")
#' round(solve(A), 3)
#' B <- createS(2000, p = 4, topology = "banded", banded.n = 2)
#' round(solve(B), 3)
#' C <- createS(2000, p = 4, topology = "clique")  # The complete graph (as m = 1)
#' round(solve(C), 3)
#' D <- createS(2000, p = 4, topology = "chain")
#' round(solve(D), 3)
#'
#' ## Generate smaple covariance matrices from block topologies:
#' C3 <- createS(2000, p = 10, topology = "clique", m = 3)
#' round(solve(C3), 1)
#' C5 <- createS(2000, p = 10, topology = "clique", m = 5)
#' round(solve(C5), 1)
#'
#' ## Can also return the precision matrix to see what happens
#' ## m = 2 blocks, each "banded" with 4 off-diagonal bands
#' round(createS(1, 12, "banded", m = 2, banded.n = 4, precision = TRUE), 2)
#'
#' ## Simulation using graph-games
#' round(createS(1, 10, "small-world", precision = TRUE), 2)
#' round(createS(1, 5, "scale-free", precision = TRUE), 2)
#' round(createS(1, 5, "random-graph", precision = TRUE), 2)
#'
#' ## Simulation using inverse Wishart distributed class covariance
#' ## Low class homogeneity
#' createS(n = c(10,10), p = 5, "banded", invwishart = TRUE, nu = 10)
#' ## Extremely high class homogeneity
#' createS(n = c(10,10), p = 5, "banded", invwishart = TRUE, nu = 1e10)
#'
#' # The precision argument can again be used to see the actual realised class
#' # precision matrices used when invwishart = TRUE.
#'
#' # The Plist argument is used to reuse old precision matrices or
#' # user-generated ones
#' P <- createS(n = 1, p = 5, "banded", precision = TRUE)
#' lapply(createS(n = c(1e5, 1e5), p = 5, Plist = list(P, P+1)), solve)
#'
#' @export
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



#' Download KEGG pathway
#'
#' Download information and graph object of a given pathway from the Kyoto
#' Encyclopedia of Genes and Genomes (KEGG) database.
#'
#' Usage of this function requires an internet connection.  The igraph objects
#' can be obtained with \code{igraph::igraph.from.graphNEL}.  The moral graph
#' can be obtained with \code{gRbase::moralize}. To obtain the adjacency matrix,
#' use \code{gRbase::as.adjMAT} or \code{igraph::get.adjacency}
#'
#' @param kegg.id A \code{character} giving the KEGG ID, e.g. \code{"map04210"},
#'   \code{"map04064"}, \code{"map04115"}. Can be prefixed with \code{"path:"}.
#'
#' @return Returns a \code{list} with entries: \item{df}{A \code{data.frame}
#'   description of the KEGG pathway.} \item{graph}{The KEGG pathway represented
#'   as a \code{graphNEL} object.}
#'
#' @note It is currently necessary to \code{require("KEGGgraph")} (or
#'   \code{require("KEGGgraph")}) due to a bug in \pkg{KEGGgraph}.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{kegg.target}}
#' @references \url{https://www.genome.jp/kegg/}
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' if (require("KEGGgraph")) {
#'   getKEGGPathway("map04064")
#' }
#' }
#' @export
getKEGGPathway <- function(kegg.id) {
  if (!requireNamespace("KEGGgraph", quietly=TRUE)) {
    stop("This function requires the bioconductor package 'KEGGgraph' and its.",
         "\ndependencies. To use it, install KEGGgraph to use it by running:\n",
         'source("http://bioconductor.org/biocLite.R")\n',
         'biocLite("KEGGgraph")\n')
  }

  # Download
  tmp.file <- paste0(tempfile(), ".kgml")
  kgml <- KEGGgraph::retrieveKGML(gsub("path:", "", kegg.id), organism="hsa",
                                  destfile = tmp.file, method = "internal")

  # Pathway data.frame information
  df        <- KEGGgraph::parseKGML2DataFrame(tmp.file)
  df$to.e   <- KEGGgraph::translateKEGGID2GeneID(df$to)
  df$from.e <- KEGGgraph::translateKEGGID2GeneID(df$from)

  # Pathway igraph and graphNEL
  graph <- KEGGgraph::KEGGpathway2Graph(KEGGgraph::parseKGML(tmp.file))

  # Make sure that the graph is simple
  # A confirmed bug in KEGGgraph and now corrected in the development branch.
  # The following lines is a temporary work-around
  graph <- igraph.to.graphNEL(simplify(igraph.from.graphNEL(graph)))

  return(list(df = df, graph = graph))
}



.parents <- function(node, graph) {
  ##############################################################################
  # - Function for extracting the parents of a node
  # - node  > A characther giving the node name in the graph
  # - graph > The graph
  #
  # NOTES:
  # - Alternative to gRbase::parents()
  ##############################################################################

  if (!requireNamespace("graph")) {
    stop("package 'graph' needed for this function.")
  }

  is.child <- sapply(graph::edges(graph), function(n) node %in% n)
  return(graph::nodes(graph)[is.child])
}



#' Construct target matrix from KEGG
#'
#' Construct a target matrix by combining topology information from the Kyoto
#' Encyclopedia of Genes and Genomes (KEGG) database and pilot data.
#'
#' The function estimates the precision matrix based on the topology given by
#' the KEGG database.  Requires a connection to the internet.
#'
#' @param Y The complete observation matrix of observations with variables in
#'   columns. The column names should be on the form e.g.  \code{"hsa:3988"}
#'   ("\code{<organism>:<Entrez id>}"). It can however also be just the Entrez
#'   id with or without the post-fixed \code{"_at"} and then the specified
#'   \code{organism} will be assumed.
#' @param kegg.id A \code{character} giving the KEGG ID, e.g. \code{"map04210"},
#'   \code{"map04064"}, or \code{"map04115"}.
#' @param method The method for estimating the non-zero entries moralized graph
#'   of the KEGG topology.  Currently, only \code{"linreg"} is implemented.
#' @param organism A \code{character} giving the organism, the default is
#'   \code{"hsa"} (homo-sapiens).
#' @param graph A \code{graphNEL} object specifying the topology of the pathway.
#'   Can be used to avoid repeatedly downloading the information.
#'
#' @return Returns a target \code{matrix} with size depending on the
#'   \code{kegg.id}.
#'
#' @note It is currently nessesary to \code{require("KEGGgraph")} (or
#'   \code{require("KEGGgraph")}) due to a bug in \pkg{KEGGgraph}.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{getKEGGPathway}}, \code{\link{default.target}}, and
#'   \code{\link{default.target.fused}}
#'
#' @references \url{https://www.genome.jp/kegg/}
#'
#' @examples
#' \dontrun{
#' if (require("KEGGgraph")) {
#' kegg.g <- getKEGGPathway("map04115")$graph
#'
#' # Create some toy data with the correct names
#' Y <- createS(n = 10, p = numNodes(kegg.g), dataset = TRUE)
#' colnames(Y) <- nodes(kegg.g)
#'
#' T <- kegg.target(Y, "map04115")
#' print(T[1:10, 1:10])
#' }
#' }
#'
#' @export kegg.target
kegg.target <- function(Y, kegg.id, method = "linreg", organism = "hsa",
                        graph = getKEGGPathway(kegg.id)$graph) {
  method <- match.arg(method)
  stopifnot(length(organism) == 1L)

  if (!requireNamespace("KEGGgraph") && !requireNamespace("graph")) {

    stop("This function requires the bioconductor package 'KEGGgraph' and its.",
         "\ndependencies. To use it, install KEGGgraph to use it by running:\n",
         'source("http://bioconductor.org/biocLite.R")\n',
         'biocLite("KEGGgraph")\n')

  }

  # Check input

  correct.format <- grepl("^([[:alpha:]]+:)?[0-9]+(_at)?$", colnames(Y))
  s <- sum(!correct.format)
  if (s > 0) {
    wmsg <- paste("Found %d colnames of Y which incorrectly formatted.",
                  "They should be on the form <organism>:<entrez id> or",
                  "<entrez id> optionally be postfixed with '_at'")
    warning(sprintf(wmsg, s))
  }

  splt <- sapply(strsplit(colnames(Y), ":"),
                 function(x) if (length(x)==2) x[1] else organism)
  if (any(splt != organism)) {
    stop("The prefix does not always match the specified organism")
  }

  #
  # Download pathway and graphNEL object
  #

  stopifnot(is(graph, "graphNEL"))
  if (!is.dag(igraph.from.graphNEL(graph))) {
    warning("The graph obtained from KEGG is not acyclic. Results are only",
            " approximate.")
  }

  # Try to correct colnames (and save the old ones)
  colnames.org <- colnames(Y)
  colnames(Y) <- gsub("_at$", "", colnames(Y))  # Remove any _at post-fix,
  colnames(Y) <- gsub(paste0("^", organism, ":"), "", colnames(Y)) # rm prefix
  colnames(Y) <- paste0(organism, ":", colnames(Y)) # Put prefix back on all

  # Determine nodes/variables both on array and in pathway and subset
  common <- intersect(graph::nodes(graph), colnames(Y))
  ind <- match(common, colnames(Y))

  if (length(common) == 0) {
    stop("There were no gene IDs in pathway and supplied data. ",
         "Check that the column names are correctly formatted.")
  }
  g <- graph::subGraph(common, graph) #=removeNode(setdiff(nodes(G),common),G)
  Ysub <- Y[, common]

  # Center the data
  Ysub <- scale(Ysub, center = TRUE, scale = FALSE)

  if (method == "linreg") {

    # Intitialize precision matrix
    p <- graph::numNodes(g)
    prec <- matrix(0, p, p)
    rownames(prec) <- colnames(prec) <- common
    for (node in graph::nodes(g)) {

      pa.node <- .parents(node, g)
      # fit <- lm(Ysub[, node] ~ Ysub[, pa.node])  # Alternative computation
      # tausq <- summary(fit)$sigma^2
      fit <- lm.fit(x = cbind(Intercept = 1, Ysub[, pa.node, drop = FALSE]),
                    y = Ysub[, node])
      tausq <- (sum(fit$residuals^2)/(nrow(Ysub) - fit$rank))
      beta <- coef(fit)

      # Update entries [node, node]
      prec[node, node] <- prec[node, node] + 1/tausq

      # Update entries [node, pa(node)]
      prec[node, pa.node] <- prec[node, pa.node] - beta[pa.node]/tausq
      prec[pa.node, node] <- prec[node, pa.node]

      # Update entries [pa(node), pa(node)]
      prec[pa.node, pa.node] <-
        prec[pa.node, pa.node] + tcrossprod(beta[pa.node])/tausq

      if (anyNA(prec)) {
        stop("NAs were introduced")
      }
    }

    if (!isSymmetricPD(prec)) {
      warning("Constructed target is not symmetric PD")
    }

    rownames(prec) <- colnames(prec) <- colnames.org[ind]
    return(prec)

  } else {

    stop("No other methods are currently implemented")

  }

}



#' Compute the pooled covariance or precision matrix estimate
#'
#' Compute the pooled covariance or precision matrix estimate from a \code{list}
#' of covariance matrices or precision matrices.
#'
#' When \code{mle} is \code{FALSE} the given covariance/precision matrices is
#' assumed to have been computed using the denominator \code{ns[i] - 1}. Hence,
#' the sum of all \code{ns} minus \eqn{G} is used a the denominator of the
#' pooled estimate. Conversely, when \code{mle} is \code{TRUE} the total sum of
#' the sample sizes \code{ns} is used as the denominator in the pooled estimate.
#'
#' The function \code{pooledP} is equivalent to a wrapper for \code{pooledS}.
#' That is, it inverts all the precision matrices in \code{Plist}, applies
#' \code{pooledS}, and inverts the resulting matrix.
#'
#' @param Slist A \code{list} of length \eqn{G} of \code{numeric} covariance
#'   matrices of the same size.
#' @param ns A \code{numeric} vector for length \eqn{G} giving the sample sizes
#'   in the corresponding entries of \code{Slist}
#' @param mle \code{logical}. If \code{TRUE}, the (biased) MLE is given. If
#'   \code{FALSE}, the biased corrected estimate is given. Default is
#'   \code{TRUE}.
#' @param subset \code{logical} vector of the same length as \code{Slist} giving
#'   the classes to pool. Default is all classes.
#' @param Plist A \code{list} of length \eqn{G} of invertible \code{numeric}
#'   precision matrices of the same size.
#'
#' @return \code{pooledS} returns the pooled covariance matrix, that is a
#'   \code{numeric} matrix with the same size as the elements of \code{Slist}.
#'   Similarly, \code{pooledP} returns the pooled precision matrix, i.e. a
#'   \code{numeric} matrix with the same size as the elements of \code{Plist}.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @examples
#' ns <- c(4, 6, 8)
#' Slist <- createS(ns, p = 6)
#'
#' pooledS(Slist, ns)
#' pooledS(Slist, ns, mle = FALSE)
#'
#' # Pool the first two classes only, leave out the remaning
#' pooledS(Slist, ns, subset = c(TRUE, TRUE, FALSE))
#' pooledS(Slist, ns, subset = ns > 5) # Pool studies with sample size > 5
#'
#' # Pooled precision matrices
#' ns <- c(7, 8, 9)
#' Plist <- lapply(createS(ns, p = 6), solve)
#' pooledS(Plist, ns)
#'
#' @export
pooledS <- function(Slist, ns, subset = rep(TRUE, length(ns)), mle = TRUE) {
  # Check inputs
  mle <- as.logical(mle)
  subset <- as.logical(subset)
  if (any(is.na(mle))) {
    stop("mle could not be coerced to a logical")
  }
  if (any(is.na(subset))) {
    stop("subset could not be coerced to a logical")
  }
  stopifnot(is.list(Slist) && length(Slist) == length(ns))
  stopifnot(is.logical(mle) && length(mle) == 1)
  stopifnot(is.logical(mle) && length(subset) == length(Slist))
  if (!any(subset)) {
    stop("argument subset must contain at least one TRUE entry.")
  }

  # Subsetting
  Slist <- Slist[subset]
  ns <- ns[subset]

  # Compute estimate
  ans <- .armaPooledS(Slist = Slist, ns = ns, mle = mle)
  dimnames(ans) <- dimnames(Slist[[1]])

  return(ans)
}



#' @rdname pooledS
#' @export
pooledP <- function(Plist, ns, subset = rep(TRUE, length(ns)), mle = TRUE) {
  # Check inputs
  mle <- as.logical(mle)
  subset <- as.logical(subset)
  if (any(is.na(mle))) {
    stop("mle could not be coerced to a logical")
  }
  if (any(is.na(subset))) {
    stop("subset could not be coerced to a logical")
  }
  stopifnot(is.list(Plist) && length(Plist) == length(ns))
  stopifnot(is.logical(mle) && length(mle) == 1)
  stopifnot(is.logical(mle) && length(subset) == length(Plist))
  if (!any(subset)) {
    stop("argument subset must contain at least one TRUE entry.")
  }

  # Subsetting
  Plist <- Plist[subset]
  ns <- ns[subset]

  # Compute estimate
  ans <- .armaPooledP(Plist = Plist, ns = ns, mle = mle)
  dimnames(ans) <- dimnames(Plist[[1]])

  return(ans)
}



#' Fused Kullback-Leibler divergence for sets of distributions
#'
#' Function calculating the Kullback-Leibler divergence between two sets of
#' multivariate normal distributions. In other words, it calculates a weigthed
#' mean of Kullback-Leibler divergences between multiple paired normal
#' distributions.
#'
#' @param MtestList A \code{list} of mean vectors of the approximating
#'   multivariate normal distribution for each class. Assumed to be zero vectors
#'   if not supplied.
#' @param MrefList A \code{list} of mean vectors of the reference multivariate
#'   normal distribution for each class. Assumed to be zero vectors if not
#'   supplied.
#' @param StestList A \code{list} of covariance matrices of the approximating
#'   multivariate normal distribtuion for each class. Usually a \code{list} of
#'   sample covariance matrices.
#' @param SrefList A \code{list} of covariance matrices of the references
#'   multivariate normal distribtuion for each class. Usually a \code{list} of
#'   the population or reference covariance matrices.
#' @param ns a \code{numeric} of the same length as the previous arguments
#'   giving the sample sizes. Used as weights in the weighted mean.
#' @param symmetric a \code{logical} indicating if original symmetric version of
#'   KL divergence should be calculated.
#'
#' @return Function returns a \code{numeric} representing the (optionally
#'   symmetric) fused Kullback-Leibler divergence.
#'
#' @author Anders Ellern Bilgrau, Wessel N. van Wieringen, Carel F.W. Peeters
#'   <cf.peeters@@vumc.nl>
#'
#' @seealso \code{\link{KLdiv}}
#'
#' @examples
#' # Create some toy data
#' n <- c(40, 60, 80)
#' p <- 10
#' Stest <- replicate(length(n), diag(p), simplify = FALSE)
#' Sref <- createS(n, p = p)
#'
#' KLdiv.fused(StestList = Stest, SrefList = Sref, ns = n, symmetric = FALSE)
#' KLdiv.fused(StestList = Stest, SrefList = Sref, ns = n, symmetric = TRUE)
#'
#' @export KLdiv.fused
KLdiv.fused <- function(MtestList, MrefList, StestList, SrefList, ns,
                        symmetric = FALSE) {

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



.LL.fused <- function(Slist, Plist, ns){
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



.PLL.fused <- function(Slist, Plist, ns, Tlist, lambda){
  ##############################################################################
  # - Function that computes the value of the (negative) penalized combined
  #   log-likelihood
  # - Slist   > A list of G sample covariance matrices for each class
  # - Plist   > A list of G precision matrices with the same size as Slist.
  #             Possibly regularized inverse of covariance matrices.
  # - ns      > A vector of sample sizes of the same length as Slist.
  # - Tlist   > A list of G p.d. target matrices
  # - lambda  > A non-negative symmetric G by G penalty matrix
  ##############################################################################

  penalty <- 0
  for (g1 in seq_along(Slist)) {
    for (g2 in seq_len(g1)) {
      if (g1 == g2) { # Ridge penalty
        penalty <- penalty +
          lambda[g1, g1]*.FrobeniusLoss(Slist[[g1]], Tlist[[g1]])
      } else {  # Fused contribution
        penalty <- penalty +
          lambda[g1, g2]*.FrobeniusLoss(Slist[[g1]] - Tlist[[g1]],
                                        Slist[[g2]] - Tlist[[g2]])
      }
    }
  }
  penalty <- penalty/2

  ans <- .LL.fused(Slist, Plist, ns) + penalty
  return(ans)
}




##------------------------------------------------------------------------------
##
## Functions for the fused ridge estimator
##
##------------------------------------------------------------------------------

.init.ridgeP.fused <- function(Slist, ns, Tlist, lambda, ...) {
  ##############################################################################
  # - Internal function for selecting initial values for Plist
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A numeric non-negative symmetric G by G penalty matrix giving
  #             the penalties of the fused ridge estimator. The diagonal entries
  #             correspond to the class ridge penalites. The off-diagonal
  #             entries, lambda[g1, g2] say, determine the retainment of
  #             similarities between estimates in classes g1 and g2.
  #             If lambda is a single number, a diagonal penalty with lambda in
  #             the diagonal is used (lambda*diag(G)).
  #             If lambda is supplied as a numeric vector of two numbers,
  #             the first is used as a common ridge penalty and the second
  #             as a common fusion penalty.
  # - ...     > Arguments passed to .armaRidgeP
  ##############################################################################

#   Spool <- pooledS(Slist, ns)
#   if (length(unique(Tlist)) == 1L) { # If all targets equal
#
#     init.Plist <-
#       .armaRidgeP(Spool, target = Tlist[[1]], .trace(lambda)/sum(ns), ...)
#     init.Plist <- replicate(length(ns), init.Plist, simplify = FALSE)
#
#   } else {
#
#     init.Plist <- vector("list", length(ns))
#     for (i in seq_along(ns)) {
#       init.Plist[[i]] <-
#         .armaRidgeP(Spool, target = Tlist[[i]], lambda[i,i]/ns[i], ...)
#     }
#
#   }

  init.Plist <- default.target.fused(Slist, ns, type = "DAIE")

  names(init.Plist) <- names(Slist)
  return(init.Plist)
}



#' Fused ridge estimation
#'
#' Performs fused ridge estimation of multiple precision matrices in cases
#' where multiple classes of data is present for given a penalty matrix.
#'
#' Performs a coordinate ascent to find the maximum likelihood of the fused
#' likelihood problem for a given ridge penalty \eqn{lambda} and fused penalty
#' matrix \eqn{Lambda_f}.
#'
#' @param Slist A \code{list} of length \eqn{G} of covariance matrices, i.e.
#' square, symmetric \code{numeric} matrices of the same size.  The \eqn{g}th
#' matrix should correspond to the \eqn{g}th class.
#' @param ns A \code{numeric} vector of sample sizes on which the matrices in
#' \code{Slist} are based.  I.e. \code{ns[g]} correspond to \code{Slist[[g]]}.
#' @param Tlist A \code{list} of length \eqn{G} of \code{numeric} p.d. target
#' matrices corresponding to the matrices in \code{Slist}.  If not supplied,
#' the default is given by \code{\link{default.target}}.
#' @param lambda The \eqn{G} by \eqn{G} penalty matrix.  That is, a symmetric,
#' non-negative \code{numeric} \code{matrix} of size \eqn{G} times \eqn{G}
#' giving the class- and pair-specific penalties.  The diagonal entries are the
#' class specific ridge penalties.  I.e. \code{lambda[i, i]} is the ridge
#' penalty for class \eqn{i}.  The off-diagonal entries are the pair-specific
#' fusion penalties.  I.e. \code{lambda[i, j]} is the fusion penalty applied on
#' the pair of classes \eqn{i} and \eqn{j}.
#'
#' Alternatively, can be supplied as a \code{numeric} of length 1 or 2.  If a
#' single number, a diagonal penalty with lambda in the diagonal is used If
#' supplied as a \code{numeric} vector of two numbers, the first is used as a
#' common ridge penalty and the second as a common fusion penalty.
#'
#' @param Plist An optional \code{list} of initial precision matrices for the
#'   fused ridge algorithm the same size as \code{Slist}.  Can be omitted.
#'   Default is the nonfused ridge precision estimate using the pooled
#'   covariance matrix corresponding to setting all fusion penalties to zero.
#' @param maxit A single \code{integer} giving the maximum number of allowed
#'   iterations.  Can be set to \code{Inf}.  If \code{maxit} is hit, a warning
#'   is given.
#' @param relative \code{logical} indicating if the convergence criterion should
#'   be on a relative scale.
#' @param verbose \code{logical}. Set to \code{TRUE} for extra output.
#' @param eps A single positive \code{numeric} giving the convergence threshold.
#'
#' @return Returns a \code{list} as \code{Slist} with precision estimates of the
#'   corresponding classes.
#'
#' @note For extreme fusion penalties in \code{lambda} the algorithm is quite
#'   sensitive to the initial values given in \code{Plist}.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{default.penalty}} \cr \code{\link{ridgeP}} for the
#'   regular ridge estimate
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#' @examples
#' # Create some (not at all high-dimensional) data on three classes
#' p <- 5  # Dimension
#' ns <- c(4, 6, 8)  # Sample sizes (K = 3 classes)
#' Slist <- createS(ns, p = p)
#' str(Slist, max.level = 2)  # The structure of Slist
#'
#' #
#' # Estimate the precisions (using the complete penalty graph)
#' #
#'
#' res1 <- ridgeP.fused(Slist, ns, lambda = c(1.3, 2.1))
#' print(res1)
#'
#' # The same using the penalty matrix (the diagnal is ignored)
#' mylambda  <- matrix(c(1.3, 2.1, 2.1,
#'                       2.1, 1.3, 2.1,
#'                       2.1, 2.1, 1.3), 3, 3, byrow = TRUE)
#' res2 <- ridgeP.fused(Slist, ns, lambda = mylambda)
#' stopifnot(all.equal(res1, res2))
#'
#'
#' #
#' # Estimate the precisions (using a non-complete penalty graph)
#' #
#'
#' # Say we only want to shrink pairs (1,2) and (2,3) and not (1,3)
#' mylambda[1,3] <- mylambda[3,1] <- 0
#' print(mylambda)
#' res3 <- ridgeP.fused(Slist, ns, lambda = mylambda)
#' # which similar to, but not the same as res1 and res2.
#'
#'
#' #
#' # Using other custom target matrices
#' #
#'
#' # Construct a custom target list
#' myTlist <- list(diag(p), matrix(1, p, p), matrix(0, p, p))
#' res4 <- ridgeP.fused(Slist, ns, Tlist = myTlist, lambda = c(1.3, 2.1))
#' print(res4)
#'
#' # Alternative, see ?default.target.fused
#' myTlist2 <- default.target.fused(Slist, ns, type = "Null")  # For the null target
#' res5 <- ridgeP.fused(Slist, ns, Tlist = myTlist2, lambda = c(1.3, 2.1))
#' print(res5)
#'
#' @export ridgeP.fused
ridgeP.fused <- function(Slist,
                         ns,
                         Tlist = default.target.fused(Slist, ns),
                         lambda,
                         Plist,
                         maxit = 100L,
                         verbose = TRUE,
                         relative = TRUE,
                         eps = sqrt(.Machine$double.eps)) {

  stopifnot(length(Slist) == length(Tlist))
  G <- length(Slist)  # Number of groups

  # Hande special imputs of lambda
  if (!is.matrix(lambda) && is.numeric(lambda)) {
    if (length(lambda) == 1) {
      lambda <- diag(lambda, G)
    } else if (length(lambda) == 2) {
      tmp <- matrix(lambda[2], G, G)
      diag(tmp) <- lambda[1]
      lambda <- tmp
    } else {
      stop("If lambda is not a numeric matrix, it should have length 1 or 2.")
    }

  }
  if (!isSymmetric(lambda, check.attributes = FALSE)) {
    stop("lambda must be symmetric")
  }
  if (!all(lambda >= 0)) {
    stop("entries of lambda must be non-negative")
  }
  if (!all(diag(lambda) > 0)) {
    stop("the diagonal of lambda must be strictly postive")
  }

  # Initialize estimates with the regular ridges from the pooled covariance
  if (missing(Plist)) {
    Plist <- .init.ridgeP.fused(Slist, ns = ns, Tlist = Tlist, lambda = lambda)
  }
  stopifnot(length(Slist) == length(Plist))

  # Overwrite the starting estimate with the fused estimate
  Plist <- .armaRidgeP.fused(Slist = Slist, ns = ns, Tlist = Tlist,
                             lambda = lambda, Plist = Plist, maxit = maxit,
                             eps = eps, relative = relative, verbose = verbose)

  # Keep dimnames and names
  for (g in seq_along(Slist)) {
    dimnames(Plist[[g]]) <- dimnames(Slist[[g]])
  }
  names(Plist) <- names(Slist)

  return(Plist)
}




##------------------------------------------------------------------------------
##
## LOOCV (and approximation) for the fused setting
##
##------------------------------------------------------------------------------

.fcvl <- function(lambda, Ylist, Tlist, init.Plist, hotstart = FALSE, ...) {
  ##############################################################################
  # - (Internal) Computes the fused leave-one-out cross-validation loss for
  #   given penalty matrix
  # - lambda     > The G by G penalty matrix.
  # - Ylist      > A list of length G of matrices of observations with samples
  #                in the rows and variables in the columns. A least 2
  #                samples (rows) are needed in each entry.
  # - Tlist      > A list of length G of target matrices the same size
  #                as those of Plist. Default is given by default.target.
  # - init.Plist > Initial estimates used in .armaRidgeP.fused
  # - hotstart   > If TRUE, the latest estimate is used for each left out
  # - ...        > Arguments passed to .armaRidgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  Slist.org <- lapply(Ylist, covML)

  # If Plist is not supplied
  if (missing(init.Plist)) {
    init.Plist <- .init.ridgeP.fused(Slist.org, ns.org, Tlist, lambda)
  }

  slh <- numeric(sum(ns.org))  # To store LOOCV losses for each sample
  j <- 1
  for (g in seq_len(G)) {
    ns <- ns.org        # "Reset" number of samples in each group
    ns[g] <- ns[g] - 1  # Update sample size in g'th group
    this.Slist <- Slist.org
    for (i in seq_len(ns.org[g])) {
      this.Slist[[g]] <-
        covML(Ylist[[g]][-i, , drop = FALSE])
        # crossprod(Ylist[[g]][-i, , drop = FALSE])/ns[g]

      this.Plist <- .armaRidgeP.fused(Slist = this.Slist, ns = ns,
                                      Tlist = Tlist, lambda = lambda,
                                      Plist = init.Plist, verbose = FALSE, ...)

      Sig <- crossprod(Ylist[[g]][i,  , drop = FALSE])
      slh[j] <- ns[g]*.LL(Sig, this.Plist[[g]])

      if (hotstart) {
        init.Plist <- this.Plist
      }

      j <- j + 1
    }
  }

  return(mean(slh))
}



.kfcvl <- function(lambda, Ylist, Tlist, init.Plist, k, ...) {
  ##############################################################################
  # - (Internal) Computes the k-fold fused cross-validation loss for a penalty
  #   matrix. The data for each class is divided into k parts. The first part
  #   in each class is left out, the fused estimate is computed based on the
  #   remaning, and the loss is computed. Then this is repeated for the remaning
  #   parts.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > Initial estimates
  # - k       > The fold of the cross validation. I.e. k is the number of
  #             roughly equally sized parts the samples for each class are
  #             partitioned into.
  # - ...     > Arguments passed to .armaRidgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  if (min(ns.org) < k) {
    stop("The least class sample size is less than the specified k = ", k)
  }

  # If Plist is not supplied
  if (missing(init.Plist)) {
    Slist.org  <- lapply(Ylist, covML)
    init.Plist <- .init.ridgeP.fused(Slist.org, ns.org, Tlist, lambda)
  }

  # Split each class into k equally sized parts
  parts <- lapply(ns.org, function(n) sample(ceiling(k*seq_len(n)/n)))

  # Run through all k splits in each class
  slh <- matrix(0, k, G)
  for (i in seq_len(k)) {
    # Pick out ALL BUT the i'th fold in each class and compute estimate
    Ylist.i <- mapply(function(x, ind) x[ind != i, , drop = FALSE],
                      Ylist, parts, SIMPLIFY = FALSE)
    ns.i    <- sapply(Ylist.i, nrow)
    Slist.i <- lapply(Ylist.i, covML)
    Plist.i <- .armaRidgeP.fused(Slist = Slist.i, ns = ns.i, Tlist = Tlist,
                                 lambda = lambda, Plist = init.Plist,
                                 verbose = FALSE, ...)

    # Evaluate estimate with left out data:
    for (g in seq_len(G)) {
      Ylist.ig <- Ylist[[g]][parts[[g]] == i, , drop = FALSE]
      nig <- nrow(Ylist.ig)
      Sig <- crossprod(Ylist.ig)/nig
      slh[i, g] <- .LL(Sig, Plist.i[[g]])
    }

  }

  return(mean(slh))
}



.sfcvl <- function(lambda, Ylist, Tlist, Plist, ...) {
  ##############################################################################
  # - (Internal) Computes the "special" LOOCV loss for given penalty parameters
  # - Only updates the class estimate in which the sample is left out.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns. A least 2
  #             samples (rows) are needed in each entry.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > Initial estimates
  # - ...     > Arguments passed to .armaRidgeP.fused and ridgeP.fused
  ##############################################################################

  G <- length(Ylist)
  ns.org <- sapply(Ylist, nrow)
  Slist.org <- lapply(Ylist, covML)

  # If Plist is not supplied
  if (missing(Plist)) {
    Plist.org <-  ridgeP.fused(Slist = Slist.org, ns = ns.org, Tlist = Tlist,
                               lambda = lambda, verbose = FALSE, ...)
  } else {
    Plist.org <- Plist
  }

  slh <- numeric(sum(ns.org))
  j <- 1
  for (g in seq_len(G)) {
    ns <- ns.org        # "Reset" number of samples in each group
    ns[g] <- ns[g] - 1  # Update sample size in g'th group

    for (i in seq_len(ns.org[g])) {
      Plist <- Plist.org
      Slist <- Slist.org
      Slist[[g]] <- covML(Ylist[[g]][-i, , drop = FALSE])

      # Update only the estimate in group "g".
      # Note these exported C++ functions are index from g = 0
      if (sum(lambda) < 1e50) {
        Plist[[g]] <- .armaFusedUpdateI(g0 = g - 1,  Plist = Plist,
                                        Slist = Slist, Tlist = Tlist, ns = ns,
                                        lambda = lambda)
      } else {
        Plist[[g]] <- .armaFusedUpdateIII(g0 = g - 1,  Plist = Plist,
                                          Slist = Slist, Tlist = Tlist, ns = ns,
                                          lambda = lambda)
      }

      Sig <- crossprod(Ylist[[g]][i,  , drop = FALSE])
      slh[j] <- .LL(Sig, Plist[[g]])
      j <- j +1
    }
  }

  return(mean(slh))
}



.afcvl <- function(lambda, Ylist, Tlist, Plist, ...) {
  ##############################################################################
  # - (Internal) Computes the approximate LOOCV loss for at given penalty
  #   parameters.
  # - lambda  > The G by G penalty matrix.
  # - Ylist   > A list of length G of matrices of observations with samples
  #             in the rows and variables in the columns.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist. Default is given by default.target.
  # - Plist   > A list of inital values for the parameter estimates
  # - ...     > Arguments passed to ridgeP.fused
  ##############################################################################

  ns <- sapply(Ylist, nrow)
  G <- length(ns)
  Slist <- lapply(Ylist, covML)
  Plist <- ridgeP.fused(Slist = Slist, Tlist = Tlist, ns = ns,
                       lambda = lambda, verbose = FALSE, ...)
  n.tot <- sum(ns)
  nll <- .LL.fused(Slist = Slist, Plist = Plist, ns)/n.tot
  p <- nrow(Slist[[1]])
  bias <- 0

  # Implementation 1
  for (g in seq_along(ns)) {
    for (i in seq_len(ns[g])) {
      Sig <- crossprod(Ylist[[g]][i, , drop = FALSE])
      fac1 <- diag(nrow(Sig)) - Sig %*% Plist[[g]]
      fac2 <- Plist[[g]] %*% (Slist[[g]] - Sig) %*% Plist[[g]]
      bias <- bias  + sum(fac1 * fac2)/(2*n.tot)
    }
  }

#   # Implementation 2  (SVANTE)
#   for (g in seq_along(ns)) {
#     lambdabar <- (lambda + sum(lambdaF[g,-g]))/ns[g]
#     P1 <- Plist[[g]]
#     P2 <- P1 %*% P1
#     P2mP1 <- P2 - Plist[[g]]
#     P4mP3 <- P2mP1 %*% P2
#     for (i in seq_len(ns[g])) {
#       yig <- Ylist[[g]][i, , drop = TRUE]
#       b1 <- (yig %*% P2mP1) %*% yig
#       b2 <- lambdabar*((yig %*% P4mP3) %*% yig)
#       bias <- bias + (b1 + b2)
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#
#   # Implementation 3 (VUCACIC)
#   for (g in seq_along(ns)) {
#     for (i in seq_len(ns[g])) {
#       Sig <- crossprod(Ylist[[g]][i, , drop = FALSE])
#       bias <- bias +
#         sum((solve(Plist[[g]]) - Sig)*Plist[[g]]*(Slist[[g]]-Sig)*Plist[[g]])
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#   # Implementation 4
#   qf <- function(x, A) return(colSums(x * (A %*% x)))
#   for (g in seq_along(ns)) {
#     ng <- ns[g]
#     t1 <- p*log((ng-1)/ng)
#     for (i in seq_len(ng)) {
#       yik <- Ylist[[g]][i, ]
#       yOy <- qf(yik, Plist[[g]])
#       oneMinusyOy <- 1 - yOy/ng
#       t1 - log(oneMinusyOy) + (yOy^2/oneMinusyOy - yOy)/ng
#     }
#   }
#   bias <- bias/(2*n.tot*(n.tot-1))
#
#
#   # Implementation 5
#   bias <- 0
#   for (g in seq_along(ns)) {
#     bias <- bias + ns[g]*sum(Plist[[g]]*(solve(Plist[[g]]) - Slist[[g]]))
#     bias <- bias + (1/rcond(Plist[[g]]))
#   }
#   bias <- bias/(2*n.tot)

  return(nll + bias)
}



.parseLambda <- function(lambda) {
  ##############################################################################
  # - A function to parse a character matrix that defines the class of penalty
  #   graphs and unique parameters for cross validation.
  # - Returns a data.frame of different indices for each level to be penalized
  #   equally.
  # - This data.frame is to be used to construct numeric matrices of penalties.
  # - lambda > A symmetric G by G character matrix defining the class of penalty
  #            matrices to cross validate over.
  #            Entries with NA, "" (the empty string), or "0" are
  #            interpreted as that that pair should omitted.
  #            Entries coercible to numeric are (in turn) interpreted as fixed
  ##############################################################################

  stopifnot(is.character(lambda))
  stopifnot(is.matrix(lambda))
  stopifnot(nrow(lambda) == ncol(lambda))

  # Handle special values
  lambda[is.na(lambda)] <- "0"
  lambda[lambda %in% ""] <- "0"

  parsedLambda <-
    data.frame(name = as.character(lambda), row = as.integer(row(lambda)),
               col = as.integer(col(lambda)), stringsAsFactors = FALSE)
  parsedLambda$val <- suppressWarnings({as.numeric(parsedLambda$name)})
  parsedLambda$fixed <- !is.na(parsedLambda$val)
  parsedLambda$variable <- !parsedLambda$fixed
  nf <- !parsedLambda$fixed
  parsedLambda$index[nf] <- as.numeric(as.factor(parsedLambda$name[nf]))

  u <- unique(subset(parsedLambda, select = -c(row, col)))

  attr(parsedLambda, "n.classes") <- nrow(lambda)
  attr(parsedLambda, "n.fixed") <- sum(u$fixed)
  attr(parsedLambda, "n.variables") <- nrow(u) - attr(parsedLambda, "n.fixed")
  return(parsedLambda)
}



.reconstructLambda <- function(lambdas, parsedLambda) {
  ##############################################################################
  # - Reconstruct the numeric penalty matrix lambda from a vector (lambdas)
  #   of penalties using the .parseLambda output.
  # - lambdas      > A numeric vector of the penalties. The length of lambdas
  #                  is the number of non-fixed entries in parsedLambda
  # - parsedLambda > A data.frame describing the penalty matrix.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  if (length(lambdas) != attributes(parsedLambda)$n.variables) {
    stop("The number of lambdas does not correspond with the number of",
         " non-fixed penalties given in parsedLambda")
  }

  G <- attributes(parsedLambda)$n.classes
  var <- parsedLambda$variable
  vals <- parsedLambda$val
  vals[var] <- lambdas[parsedLambda$index[var]]
  lambda <- matrix(vals, G, G)
  return(lambda)
}



.lambdasFromMatrix <- function(lambda.init, parsedLambda) {
  ##############################################################################
  # - Create the "lambdas" vector used in the optimizers from a numeric matrix.
  # - lambda.init  > A numeric matrix of the initial penalty matrix.
  # - parsedLambda > A data.frame describing the penalty matrix.
  #                  Should be the output from .parseLambda.
  ##############################################################################

  u <- parsedLambda[!duplicated(parsedLambda$index), ]
  u <- u[!u$fixed, ]
  lambdas <- numeric(attributes(parsedLambda)$n.variables)
  stopifnot(length(lambdas) == nrow(u))
  lambdas[u$index] <- lambda.init[as.matrix(subset(u, select = c(row, col)))]

  # Try to reconstruct the given matrix
  relambda.init <- .reconstructLambda(lambdas, parsedLambda)
  if (!all(relambda.init == lambda.init)) {
    warning("The fixed penalties do not agree with the specified initial ",
            "penalty matrix.")
  }
  return(lambdas)
}


#' @rdname optPenalty.fused
#' @export
optPenalty.fused.grid <-
  function(Ylist, Tlist,
           lambdas = 10^seq(-5, 5, length.out = 15),
           lambdaFs = lambdas,
           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
           k = 10,
           verbose = TRUE,
           ...) {
  ##############################################################################
  # - Cross validation for the fused ridge estimator on a grid to determine
  #   optimal lambda and lambdaF.
  # - Ylist       > A list of length G of matrices of observations with samples
  #                 in the rows and variables in the columns.
  # - Tlist       > A list of length G of target matrices the same size
  #                 as those of Plist. Default is given by default.target.
  # - lambdas     > A vector of ridge penalties
  # - lambdaFs    > A vector of fusion penalties
  # - cv.method   > The LOOCV type to use. Allowed values are LOOCV, aLOOCV,
  #                 sLOOCV, kCV for leave-one-out cross validation (LOOCV),
  #                 appproximate LOOCV, special LOOCV, and k-fold CV, resp.
  # - k           > Number of parts in k-fold CV. Only use if method is "kCV".
  # - ...         > Arguments passed to ridgeP.fused
  # - verbose     > logical. Print extra information. Defaults is TRUE.
  #
  # NOTES:
  # - The complete penalty graph is assumed (i.e. all ridge penalties
  #   equal and all fusion penalties equal)
  ##############################################################################

  cv.method <- match.arg(cv.method)

  if (missing(Tlist)) {  # If Tlist is not provided
    Tlist <- lapply(Ylist, function(Y) default.target(covML(Y)))
  }

  G <- length(Ylist)

  stopifnot(all(lambdas > 0))
  stopifnot(all(lambdaFs >= 0))
  stopifnot(all(is.finite(lambdas)))
  stopifnot(all(is.finite(lambdaFs)))

  if (cv.method == "LOOCV") {
    cvfunc <- .fcvl
  } else if (cv.method == "aLOOCV") {
    cvfunc <- .afcvl
  } else if (cv.method == "sLOOCV") {
    cvfunc <- .sfcvl
  } else if (cv.method == "kCV") {
    cvfunc <-  function(lambda, Ylist = Ylist, Tlist = Tlist, ...) {
      .kfcvl(lambda, Ylist = Ylist, Tlist = Tlist, k = k, ...)
    }
  } else {
    stop("cv.method not implmented.")
  }

  # Calculate CV scores
  if (verbose) {
    cat("Calculating cross-validated negative log-likelihoods...\n")
  }

  slh <- matrix(NA, length(lambdas), length(lambdaFs))
  dimnames(slh) <- list("lambdas" = lambdas, "lambdaFs" = lambdaFs)
  for (l1 in seq_along(lambdas)) {
    for (l2 in seq_along(lambdaFs)) {
      # Create penalty matrix
      lambda <- matrix(lambdaFs[l2], G, G)
      diag(lambda) <- lambdas[l1]

      # Evaluate loss
      slh[l1, l2] <- cvfunc(lambda = lambda, Ylist = Ylist, Tlist = Tlist, ...)

      if (verbose) {
        cat(sprintf("lambda = %.3e (%d), lambdaF = %.3e (%d), fcvl = %.3e\n",
                   lambdas[l1],  l1, lambdaFs[l2], l2, slh[l1, l2]))
      }
    }
  }

  output <- list(lambda = lambdas, lambdaF = lambdaFs, fcvl = slh)
  class(output) <- "optPenaltyFusedGrid"
  return(output)
}



#' Print and plot functions for fused grid-based cross-validation
#'
#' Print and plot functions for the output from
#' \code{\link{optPenalty.fused.grid}} which performs a grid based
#' cross-validation (CV) search to find optimal penalty parameters. Currently,
#' only the complete penalty graph is supported.
#'
#' @param x A \code{optPenaltyFusedGrid}-object print or plot.  Usually the
#'   output of \cr \code{\link{optPenalty.fused.grid}}.
#' @param add.text A \code{logical} value controlling if the text should be
#'   added to the plot.
#' @param add.contour A \code{logical} value controlling if the contour lines
#'   should be added to the plot.
#' @param col A \code{character} vector of colours used in the image plot.
#' @param \dots Arguments passed on.  In \code{print.optPenaltyFusedGrid} the
#'   arguments are passed to \code{print.matrix}.  In
#'   \code{plot.optPenaltyFusedGrid} are passed to the standard \code{plot}
#'   function.
#'
#' @return Invisibly returns the object (\code{x}).
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{optPenalty.fused.grid}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#' @export
print.optPenaltyFusedGrid <- function(x, ...) {
  with(x, print(fcvl))
  return(invisible(x))
}


#' @rdname print.optPenaltyFusedGrid
#' @export
plot.optPenaltyFusedGrid <- function(x, add.text = TRUE, add.contour = TRUE,
                                     col = rainbow(100, end = 0.8), ...) {
  with(x, {
    image(log(lambda), log(lambdaF), fcvl, col = col)
    if (add.contour) {
      contour(log(lambda), log(lambdaF), log(fcvl - min(fcvl) + 1), add = TRUE,
              nlevels = 15, col = "White", drawlabels = FALSE)
    }
    cols <- with(x, ifelse(fcvl == min(fcvl), "red", "black"))
    if (add.text) {
      text(log(lambda)[c(row(fcvl))], log(lambdaF)[c(col(fcvl))],
           sprintf("%0.1f", fcvl), cex = 0.7, col = cols)
    }
  })
  return(invisible(x))
}

#' @rdname optPenalty.fused
#' @importFrom stats optim
#' @export
optPenalty.fused.auto <-
  function(Ylist,
           Tlist,
           lambda,
           cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
           k = 10,
           verbose = TRUE,
           lambda.init,
           maxit.ridgeP.fused = 1000,
           optimizer = "optim",
           maxit.optimizer = 1000,
           debug = FALSE,
           optim.control = list(trace = verbose, maxit = maxit.optimizer),
           ...) {

  cv.method <- match.arg(cv.method)
  G <- length(Ylist)

  if (missing(lambda)) {  # Handle missing lambda
    lambda <- matrix("fusion", G, G)
    diag(lambda) <- "ridge"
  }

  # Interpret given lambda
  parsedLambda <- .parseLambda(lambda)

  if (verbose) {  # Report interpretation
    n.fixed     <- attributes(parsedLambda)$n.fixed
    n.variables <- attributes(parsedLambda)$n.variables
    nonfix <- with(parsedLambda, paste(unique(name[!fixed]), collapse = ", "))
    fixed  <- with(parsedLambda, paste(unique(name[ fixed]), collapse = ", "))
    message("Found ", n.fixed + n.variables, " unique penalties of which ",
            n.fixed, " are interpreted as fixed and ", n.variables,
            " are to be determined by ", cv.method, ".\n",
            "Non-fixed parameters: ", nonfix,
            "\nFixed parameters: ", ifelse(n.fixed, fixed, "<none>"))
  }

  # Determine what loss function to use
  if (cv.method == "LOOCV") {
    cvfunc <- .fcvl
  } else if (cv.method == "aLOOCV") {
    cvfunc <- .afcvl
  } else if (cv.method == "sLOOCV") {
    cvfunc <- .sfcvl
  } else if (cv.method == "kCV") {
    cvfunc <-  function(lambda, Ylist = Ylist, Tlist = Tlist, ...) {
      .kfcvl(lambda, Ylist = Ylist, Tlist = Tlist, k = k, ...)
    }
  } else {
    stop("cv.method not implmented.")
  }

  # Construct optim/nlm objective function, parameters are assumed on log-scale
  cvl <- function(lambdas, ...) {
    elambdas <- exp(lambdas)
    lambda <- .reconstructLambda(elambdas, parsedLambda)
    return(cvfunc(lambda = lambda, Ylist = Ylist, Tlist = Tlist,
                  maxit = maxit.ridgeP.fused, ...))
  }

  if (missing(lambda.init)) {
    # Get somewhat sensible starting value for non-fixed diagonal entries
    # (ridge penalties) and by choosing off-diag lambda to be zero.
    f <- function(x) {
      lambdas <- suppressWarnings({.lambdasFromMatrix(diag(exp(x), G),
                                                      parsedLambda)})
      return(cvl(lambdas))
    }
    st <- optimize(f, lower = -20, upper = 20)$minimum
    lambdas.st <- suppressWarnings({log(.lambdasFromMatrix(diag(exp(st), G),
                                                           parsedLambda))})
    lambdas.st[lambdas.st == -Inf] <- -40
    lambdas.st[lambdas.st ==  Inf] <-  40

  } else {

    if (is.matrix(lambda.init) && is.numeric(lambda.init) &&
        isSymmetric(lambda.init)) {
      lambdas.st <- log(.lambdasFromMatrix(lambda.init, parsedLambda))
      lambdas.st[lambdas.st == -Inf] <- -40
      lambdas.st[lambdas.st ==  Inf] <-  40

    } else {
      stop("The supplied inital parameters must be a symmetric numeric matrix")
    }
  }

  if (optimizer == "optim") {

    ans <- optim(lambdas.st, fn = cvl, ..., control = optim.control)
    par <- ans$par
    val <- ans$value

    args <- list(...)
    if (!is.null(args$method)) {
      if (args$method != "SANN" && ans$convergence == 1) {
        warning("Iteration limit of optim had been reached.")
      }
      if (args$method == "Nelder-Mead" && ans$convergence == 10) {
        warning("Degeneracy of the Nelder-Mead simplex.")
      }
    }

  } else if (optimizer == "nlm") {

    ans <- nlm(cvl, lambdas.st, iterlim = maxit.optimizer, ...)
    par <- ans$estimate
    val <- ans$minimum

  }

  # Format optimal values
  opt.lambdas <- exp(par)
  opt.lambda <- .reconstructLambda(opt.lambdas, parsedLambda)
  dimnames(opt.lambda) <- dimnames(lambda)

  # Construct output
  res <- list(Plist = NA,
              lambda = opt.lambda,
              lambda.unique = NA,
              value = val)

  # Compute estimate at optimal values
  res$Plist <- ridgeP.fused(Slist = lapply(Ylist, covML),
                            ns = sapply(Ylist, nrow),
                            Tlist = Tlist, lambda = res$lambda,
                            maxit = maxit.ridgeP.fused, verbose = FALSE)

  res$lambda.unique <- unique(opt.lambda[lower.tri(opt.lambda, diag = TRUE)])
  names(res$lambda.unique) <- lambda[match(res$lambda.unique, opt.lambda)]

  if (debug) {
    attr(res, "optim.debug") <- ans
  }

  return(res)
}



#' Identify optimal ridge and fused ridge penalties
#'
#' Functions to find the optimal ridge and fusion penalty parameters via
#' leave-one-out cross validation. The functions support leave-one-out
#' cross-validation (LOOCV), \eqn{k}-fold CV, and two forms of approximate
#' LOOCV. Depending on the used function, general numerical optimization or a
#' grid-based search is used.
#'
#' \code{optPenalty.fused.auto} serves a utilizes \code{\link{optim}} for
#' identifying the optimal fused parameters and works for general classes of
#' penalty graphs.
#'
#' \code{optPenalty.fused.grid} gives a grid-based evaluation of the
#' (approximate) LOOCV loss.
#'
#' @param Ylist A \code{list} of \eqn{G} matrices of data with \eqn{n_g} samples
#'   in the rows and \eqn{p} variables in the columns corresponding to \eqn{G}
#'   classes of data.
#' @param Tlist A \code{list} of \eqn{G} of p.d. class target matrices of size
#'   \eqn{p} times \eqn{p}.
#' @param lambda A symmetric \code{character} \code{matrix} encoding the class
#'   of penalty matrices to cross-validate over.  The diagonal elements
#'   correspond to the class-specific ridge penalties whereas the off-diagonal
#'   elements correspond to the fusion penalties.  The unique elements of lambda
#'   specify the penalties to determine by the method specified by
#'   \code{cv.method}.  The penalties can be fixed if they are coercible to
#'   numeric values, such as e.g. \code{"0"}, \code{"2.71"} or \code{"3.14"}.
#'   Fusion between pairs can be "left out"" using either of \code{""},
#'   \code{NA}, \code{"NA"}, or \code{"0"}.  See \code{\link{default.penalty}}
#'   for help on the construction hereof and more details.  Unused and can be
#'   omitted if \code{grid == TRUE}.
#' @param cv.method \code{character} giving the cross-validation (CV) to use.
#'   The allowed values are \code{"LOOCV"}, \code{"aLOOCV"}, \code{"sLOOCV"},
#'   \code{"kCV"} for leave-one-out cross validation (LOOCV), appproximate
#'   LOOCV, special LOOCV, and k-fold CV, respectively.
#' @param k \code{integer} giving the number of approximately equally sized
#'   parts each class is partioned into for \eqn{k}-fold CV.  Only use if
#'   \code{cv.method} is \code{"kCV"}.
#' @param verbose \code{logical}. If \code{TRUE}, progress information is
#'   printed to the console.
#' @param lambda.init A \code{numeric} penalty \code{matrix} of initial values
#'   passed to the optimizer. If omitted, the function selects a starting values
#'   using a common ridge penaltiy (determined by 1D optimization) and all
#'   fusion penalties to zero.
#' @param maxit.ridgeP.fused A \code{integer} giving the maximum number of
#'   iterations allowed for each fused ridge fit.
#' @param optimizer \code{character}. Either \code{"optim"} or \code{"nlm"}
#'   determining which optimizer to use.
#' @param maxit.optimizer A \code{integer} giving the maximum number of
#'   iterations allowed in the optimization procedure.
#' @param debug \code{logical}. If \code{TRUE} additional output from the
#'   optimizer is appended to the output as an attribute.
#' @param lambdas A \code{numeric} vector of positive ridge penalties.
#' @param lambdaFs A \code{numeric} vector of non-negative fusion penalties.
#' @param grid \code{logical.} Should a grid based search be used? Default is
#'   \code{FALSE}.
#' @param optim.control A \code{list} of control arguments for
#'   \code{\link{optim}}.
#' @param \dots For \code{optPenalty.fused}, arguments are passed to
#'   \code{optPenalty.fused.grid} or \code{optPenalty.fused.auto} depending on
#'   the value of \code{grid}.  In \code{optPenalty.fused.grid}, arguments are
#'   passed to \code{ridgeP.fused}.  In \code{optPenalty.fused.auto}, arguments
#'   are passed to the optimizer.
#'
#' @return \code{optPenalty.fused.auto} returns a \code{list}:\cr \item{Plist}{A
#'   \code{list} of the precision estimates for the optimal parameters.}
#'   \item{lambda}{The estimated optimal fused penalty matrix.}
#'   \item{lambda.unique}{The unique entries of the \code{lambda}.  A more
#'   concise overview of \code{lambda}} \item{value}{The value of the loss
#'   function in the estimated optimum.}
#'
#'   \code{optPenalty.fused.LOOCV} returns a \code{list}:\cr \item{ridge}{A
#'   \code{numeric} vector of grid values for the ridge penalty}
#'   \item{fusion}{The \code{numeric} vector of grid values for the fusion
#'   penalty} \item{fcvl}{The \code{numeric} \code{matrix} of evaluations of the
#'   loss function}
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso See also \code{\link{default.penalty}}, \code{optPenalty.LOOCV}.
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#' @examples
#' \dontrun{
#' # Generate some (not so) high-dimensional data witn (not so) many samples
#' ns <- c(4, 5, 6)
#' Ylist <- createS(n = ns, p = 6, dataset = TRUE)
#' Slist <- lapply(Ylist, covML)
#' Tlist <- default.target.fused(Slist, ns, type = "DIAES")
#'
#'
#' # Grid-based
#' lambdas <- 10^seq(-5, 3, length.out = 7)
#' a <- optPenalty.fused.grid(Ylist, Tlist,
#'                            lambdas = lambdas,
#'                            cv.method = "LOOCV", maxit = 1000)
#' b <- optPenalty.fused.grid(Ylist, Tlist,
#'                            lambdas = lambdas,
#'                            cv.method = "aLOOCV", maxit = 1000)
#' c <- optPenalty.fused.grid(Ylist, Tlist,
#'                            lambdas = lambdas,
#'                            cv.method = "sLOOCV", maxit = 1000)
#' d <- optPenalty.fused.grid(Ylist, Tlist,
#'                            lambdas = lambdas,
#'                            cv.method = "kCV", k = 2, maxit = 1000)
#'
#' # Numerical optimization (uses the default "optim" optimizer with method "BFGS")
#' aa <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "LOOCV", method = "BFGS")
#' print(aa)
#' bb <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "aLOOCV", method = "BFGS")
#' print(bb)
#' cc <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "sLOOCV", method = "BFGS")
#' print(cc)
#' dd <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "kCV", k=3, method="BFGS")
#' print(dd)
#'
#'
#' #
#' # Plot the results
#' #
#'
#' # LOOCV
#' # Get minimums and plot
#' amin  <- log(expand.grid(a$lambda, a$lambdaF))[which.min(a$fcvl), ]
#' aamin <- c(log(aa$lambda[1,1]), log(aa$lambda[1,2]))
#'
#' # Plot
#' filled.contour(log(a$lambda), log(a$lambdaF), log(a$fcvl), color = heat.colors,
#'                plot.axes = {points(amin[1], amin[2], pch = 16);
#'                             points(aamin[1], aamin[2], pch = 16, col = "purple");
#'                             axis(1); axis(2)},
#'                xlab = "lambda", ylab = "lambdaF", main = "LOOCV")
#'
#' # Approximate LOOCV
#' # Get minimums and plot
#' bmin <- log(expand.grid(b$lambda, b$lambdaF))[which.min(b$fcvl), ]
#' bbmin <- c(log(bb$lambda[1,1]), log(unique(bb$lambda[1,2])))
#'
#' filled.contour(log(b$lambda), log(b$lambdaF), log(b$fcvl), color = heat.colors,
#'                plot.axes = {points(bmin[1], bmin[2], pch = 16);
#'                             points(bbmin[1], bbmin[2], pch = 16, col ="purple");
#'                             axis(1); axis(2)},
#'                xlab = "lambda", ylab = "lambdaF", main = "Approximate LOOCV")
#'
#'
#' #
#' # Arbitrary penalty graphs
#' #
#'
#' # Generate some new high-dimensional data and a 2 by 2 factorial design
#' ns <- c(6, 5, 3, 2)
#' df <- expand.grid(Factor1 = LETTERS[1:2], Factor2 = letters[3:4])
#' Ylist <- createS(n = ns, p = 4, dataset = TRUE)
#' Tlist <- lapply(lapply(Ylist, covML), default.target, type = "Null")
#'
#' # Construct penalty matrix
#' lambda <- default.penalty(df, type = "CartesianUnequal")
#'
#' # Find optimal parameters,
#' # Using optim with method "Nelder-Mead" with "special" LOOCV
#' ans1 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
#'                          cv.method = "sLOOCV", verbose = FALSE)
#' print(ans1$lambda.unique)
#'
#' # By approximate LOOCV using optim with method "BFGS"
#' ans2 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
#'                          cv.method = "aLOOCV", verbose = FALSE,
#'                          method = "BFGS")
#' print(ans2$lambda.unique)
#'
#' # By LOOCV using nlm
#' lambda.init <- matrix(1, 4, 4)
#' lambda.init[cbind(1:4,4:1)] <- 0
#' ans3 <- optPenalty.fused(Ylist, Tlist, lambda = lambda,
#'                          lambda.init = lambda.init,
#'                          cv.method = "LOOCV", verbose = FALSE,
#'                          optimizer = "nlm")
#' print(ans3$lambda.unique)
#'
#' # Quite different results!
#'
#'
#' #
#' # Arbitrary penalty graphs with fixed penalties!
#' #
#'
#' # Generate some new high-dimensional data and a 2 by 2 factorial design
#' ns <- c(6, 5, 5, 5)
#' df <- expand.grid(DS = LETTERS[1:2], ER = letters[3:4])
#' Ylist <- createS(n = ns, p = 4, dataset = TRUE)
#' Tlist <- lapply(lapply(Ylist, covML), default.target, type = "Null")
#'
#' lambda <- default.penalty(df, type = "Tensor")
#' print(lambda)  # Say we want to penalize the pair (1,2) with strength 2.1;
#' lambda[2,1] <- lambda[1,2] <- 2.1
#' print(lambda)
#'
#' # Specifiying starting values is also possible:
#' init <- diag(length(ns))
#' init[2,1] <- init[1,2] <- 2.1
#'
#' res <- optPenalty.fused(Ylist, Tlist, lambda = lambda, lambda.init = init,
#'                         cv.method = "aLOOCV", optimizer = "nlm")
#' print(res)
#' }
#'
#' @export optPenalty.fused
optPenalty.fused <- function(Ylist, Tlist, lambda = default.penalty(Ylist),
                             cv.method = c("LOOCV", "aLOOCV", "sLOOCV", "kCV"),
                             k = 10, grid = FALSE, ...) {
  cv.method <- match.arg(cv.method)

  if (grid) {
    res <- optPenalty.fused.grid(Ylist = Ylist, Tlist = Tlist,
                                 cv.method = cv.method, k = k, ... )
  } else {
    res <- optPenalty.fused.auto(Ylist = Ylist, Tlist = Tlist,
                                 lambda = lambda,
                                 cv.method = cv.method, k = k,...)
  }
  return(res)
}




##------------------------------------------------------------------------------
##
## Automatic penalty matrix constructor
##
##------------------------------------------------------------------------------

.charAdjMat <- function(fac, name = "X", ordered = is.ordered(fac)) {
  ##############################################################################
  # - Create a complete character adjacency matrix from a factor. This function
  #   is used in the constructing the character penalty matrix in
  #   default.penalty.
  # - fac     > A factor of some length. Can be ordered.
  # - name    > A character giving the text which should appear in the adjacent
  #             entries. If not a character, the object name of fac is used.
  # - ordered > logical specifiying if fac should be interpreted as ordered.
  #
  # Examples:
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]))
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = "Y")
  #  rags2ridges:::.charAdjMat(factor(LETTERS[1:3]), name = NULL)
  #  rags2ridges:::.charAdjMat(ordered(factor(LETTERS[1:5])))
  ##############################################################################

  G <- nlevels(fac)
  if (is.character(name)) {
    lab <- name
  } else {
    lab <- deparse(substitute(fac))
  }
  if (ordered) {
    M <- matrix("", G, G)
    M[row(M) == (col(M) + 1)] <- lab
    M[row(M) == (col(M) - 1)] <- lab
  } else {
    M <- matrix(lab, G, G)
    diag(M) <- ""
  }
  rownames(M) <- colnames(M) <- levels(fac)
  return(M)
}



.char2num <- function(X) {
  ##############################################################################
  # - Create a character adjacency matrix to a numeric one
  # - X  > A character matrix where "" signify non-adjacency.
  #
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
  # - Construct the Cartesian product graph from two "character" matrices.
  # - A > A character matrix where "" signify non-adjacency.
  # - B > A character matrix where "" signify non-adjacency.
  #
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
  # - Construct the Tensor (or categorical) product graph from two "character"
  #   matrices.
  # - A > A character matrix where "" signify non-adjacency.
  # - B > A character matrix where "" signify non-adjacency.
  #
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



#' Construct commonly used penalty matrices
#'
#' Function that constructs default or commonly use penalty matrices according
#' to a (factorial) study design.  The constructed penalty matrix can be used
#' directly in \code{\link{optPenalty.fused.auto}} or serve as basis for
#' modification.
#'
#' The \code{type} gives a number of common choices for the penalty matrix:
#' \itemize{ \item \code{'Complete'} is the complete penalty graph with equal
#' penalties.  \item \code{'CartesianEqual'} corresponds to a penalizing along
#' each "direction" of factors with a common penalty. The choice is named
#' Cartesian as it is the Cartesian graph product of the complete penalty
#' graphs for the individual factors.  \item \code{'CartesianUnequal'}
#' corresponds to a penalizing each direction of factors with individual
#' penalties.  \item \code{'TensorProd'} correspond to penalizing the
#' "diagonals" only.  It is equivalent to the graph tensor products of the
#' complete graphs for each individual factor.  }
#'
#' @param G A \code{numeric} giving the number of classes. Can also be a
#'   \code{list} of length \code{G} such as the usual argument \code{Slist} from
#'   other \pkg{rags2ridges} functions.  Can be omitted if \code{df} is given.
#' @param df A \code{data.frame} with \code{G} rows and factors in the columns.
#'   Note, the columns has to be of type \code{factor}.  Can be omitted when
#'   \code{G} is given and \code{type == "Complete"}.  The factors can be
#'   ordered.
#' @param type A character giving the type of fused penalty graph to construct.
#'   Should be \code{'Complete'} (default), \code{'CartesianEqual'}, or
#'   \code{'CartesianUnequal'} or \code{'TensorProd'} or an unique abbreviation
#'   hereof. See details.
#'
#' @return Returns a \code{G} by \code{G} character matrix which specify the
#'   class of penalty graphs to be used.  The output is suitable as input for
#'   the penalty matrix used in \code{\link{optPenalty.fused.auto}}.
#'
#' @author Anders E. Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel
#'   N. van Wieringen
#'
#' @seealso \code{\link{ridgeP.fused}}, \code{\link{optPenalty.fused}},
#'   \code{\link{default.target}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#' @examples
#'   # Handling one-way designs
#'   default.penalty(2)
#'   default.penalty(4)
#'   Slist <- vector("list", 6)
#'   default.penalty(Slist)   # The function uses only the length of the list
#'   df0 <- expand.grid(Factor = c("lvl1", "lvl2"))
#'   default.penalty(df0)
#'
#'   # A more elaborate example
#'   df1 <- expand.grid(DS = c("DS1", "DS2", "DS3"), ER = c("ER+", "ER-"))
#'
#'   # Usage (various interface demonstrations)
#'   default.penalty(6, df1, type = "Complete")
#'   default.penalty(6, type = "CartesianEqual")  # GIVES WARNING
#'   default.penalty(6, df1, type = "CartesianEqual")
#'   default.penalty(Slist, df1, type = "CartesianEqual")
#'   default.penalty(6, df1, type = "CartesianUnequal")
#'   default.penalty(df1)
#'
#'   # A 2 by 2 by 2 design
#'   df2 <- expand.grid(A = c("A1", "A2"), B = c("B1", "B2"), C = c("C1", "C3"))
#'   default.penalty(df2)
#'   default.penalty(df2, type = "CartesianEqual")
#'   default.penalty(df2, type = "CartesianUnequal")
#'
#' @export default.penalty
default.penalty <- function(G, df,
                            type = c("Complete", "CartesianEqual",
                                     "CartesianUnequal", "TensorProd")) {

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

  # Make sure the levels of the factors are ordered correctly
  for (i in seq_len(ncol(df))) {
    df[[i]] <- factor(df[[i]], levels = unique(df[[i]]))
  }

  # Construct penalty matrix class
  if (type == "Complete") {

    M <- matrix("fusion", G, G)
    rownames(M) <- colnames(M) <- Reduce(":", df)
    diag(M) <- "ridge"
    return(M)

  } else if (type == "CartesianEqual" || type == "CartesianUnequal") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.cartesianProd, adj.mats)

    if (type == "CartesianEqual") {
      M[M != ""] <- "fusion"
    }
    diag(M) <- "ridge"
    return(M)

  } else if (type == "TensorProd") {

    adj.mats <- lapply(seq_along(df),
                       function(i) .charAdjMat(df[[i]], name = names(df)[i]))
    M <- Reduce(.tensorProd, adj.mats)
    diag(M) <- "ridge"
    return(M)

  } else {

    stop("type =", type, "not implemented yet!")

  }
}




##------------------------------------------------------------------------------
##
## To fuse or not to fuse --- Test H0: Omega_1 = ... = Omega_G
##
##------------------------------------------------------------------------------

.scoreStatistic <- function(Plist, Slist, ns) {
  ##############################################################################
  # - Function for computing the score statistic
  # - Plist > A list of precision matrices
  # - Slist > A list of sample covariance matrices
  # - ns    > A vector with the same length as Plist and Slist of sample sizes
  ##############################################################################

  stopifnot(length(Plist) == length(Slist))
  stopifnot(length(ns) == length(Plist))

  U <- 0
  for (g in seq_along(ns)) {
    X <- ns[g]*(Plist[[g]] - Slist[[g]])
    diag(X) <- 0.5*diag(X)
    U <- U + sum(X * (Plist[[g]] %*% X %*% Plist[[g]]))
  }
  return(U)
}



.scambleYlist <- function(Ylist) {
  ##############################################################################
  # - Function for permuting the class labels of Ylist, equivalent to
  #   scrambling/permuting all obervations.
  # - Ylist > A list of observations matrices for each class
  ##############################################################################

  ns <- sapply(Ylist, nrow)
  cl <- factor(rep(names(Ylist), ns))
  Y  <- as.data.frame(do.call(rbind, Ylist))
  out <- split(Y, sample(cl))
  out <- lapply(out, as.matrix)
  return(out)
}



#' Test the necessity of fusion
#'
#' Function for testing the null hypothesis that all population precision
#' matrices are equal and thus the necessity for the fusion penalty. Note, the
#' test performed is conditional on the supplied penalties and targets.
#'
#' The function computes the observed score statistic \eqn{U_obs} using the
#' fused ridge estimator on the given data. Next, the score statistic is
#' computed a number of times (given by \code{n.permutations}) under the
#' null-hypothesis by effectively permuting the class labels of the data.
#'
#' @param Ylist A \code{list} of length \eqn{G} of observations matrices for
#'   each class.  Variables are assumed to correspond to the columns.
#' @param Tlist A \code{list} of target matrices for each class. Should be same
#'   length as \code{Ylist}-
#' @param lambda A non-negative, symmetric \eqn{G} by \eqn{G} \code{matrix}
#'   giving the ridge and fusion penalties.
#' @param n.permutations The number of permutations to approximate the null
#'   distribution.  Default is 100. Should be increased if sufficient computing
#'   power is available.
#' @param verbose Print out extra progress information
#' @param \dots Arguments passed to \code{\link{ridgeP.fused}}.
#'
#' @return Returns a \code{list} values containing the observed test statistic
#'   and the test statistic under the null distribution.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel, N. van Wieringen
#'
#' @seealso \code{\link{ridgeP.fused}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#' @examples
#' ns <- c(10, 5, 23)
#' Ylist <- createS(ns, p = 15, topology = "banded", dataset = TRUE)
#'
#' # Use the identity target matrix for each class
#' Tlist <- replicate(length(ns), diag(15), simplify = FALSE)
#'
#' # Do the test
#' lm <- matrix(10, 3, 3)
#' diag(lm) <- 1
#' ft <- fused.test(Ylist, Tlist, lambda = lm,
#'                  n.permutations = 500)
#' print(ft)
#'
#' # Summary spits out a bit more information
#' summary(ft)
#'
#' # The returned object can alo be plotted via
#' hist(ft)
#' # or via the alias
#' plot(ft)
#'
#' # Customization and parameters work a usual:
#' hist(ft, col = "steelblue", main = "Null distribution", add.extra = FALSE,
#'      xlab = "Score statistic", freq = FALSE)
#'
#' @export
fused.test <- function(Ylist, Tlist, lambda,
                       n.permutations = 100, verbose = FALSE, ...) {
  stopifnot(length(Ylist) == length(Tlist))
  stopifnot(nrow(lambda) == length(Ylist))
  stopifnot(ncol(lambda) == length(Ylist))

  G <- length(Ylist)
  ns <- sapply(Ylist, nrow)
  n.tot <- sum(ns)
  lambda.null <- G*lambda/sum(ns)

  # Compute observed statistic
  if (verbose) {message("Computing the observed score statistic... ")}
  Slist <- lapply(Ylist, covML)
  Spool.obs <- pooledS(Slist, ns)
  Plist.obs <- list()
  for (i in seq_len(G)) {
    Plist.obs[[i]] <- .armaRidgeP(Spool.obs, target = Tlist[[i]],
                                  lambda = lambda.null[i,i])
  }

  Uobs <- .scoreStatistic(Plist = Plist.obs, Slist = Slist, ns = ns)

  # Approximate null distribution by permutation
  if (verbose) {message("Computing the score statistics under permutation... ")}
  Unull <- numeric()
  for (j in seq_len(n.permutations)) {
    Ylist.tmp <- .scambleYlist(Ylist)  # Permute class labels
    Spool.tmp <- pooledS(lapply(Ylist.tmp, covML), ns)
    Plist.null <- list()
    for (i in seq_len(G)) {
      Plist.null[[i]] <- .armaRidgeP(Spool.tmp,
                                     target = Tlist[[i]],
                                     lambda = lambda.null[i,i])
    }
    Unull[j] <- .scoreStatistic(Plist = Plist.null,
                                Slist = lapply(Ylist.tmp, covML), ns = ns)

    if (verbose && j %% 10 == 0) {
      cat(sprintf("%d of %d done\n", j, n.permutations))
    }
  }

  # Return results
  ans <- list(observed = Uobs, null.dist = Unull)
  class(ans) <- "ptest"
  return(ans)
}



#' Print and summarize fusion test
#'
#' Print and summary functions for the fusion test performed by
#' \code{\link{fused.test}}.
#'
#' @param x,object The object to print or summarize. Usually the output of
#'   \code{\link{fused.test}}.
#' @param digits An \code{integer} controlling the number of printed digits.
#' @param \dots Arguments passed on.  In \code{summary.ptest} the arguments are
#'   passed to \code{print.ptest}.  In \code{print.ptest} are passed to the
#'   standard \code{summary} function.
#'
#' @return Invisibly returns the object.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{fused.test}}, \code{\link{hist.ptest}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#' @examples
#' ns <- c(10, 5, 23)
#' Ylist <- createS(ns, p = 15, topology = "banded", dataset = TRUE)
#'
#' # Use the identity target matrix for each class
#' Tlist <- replicate(length(ns), diag(15), simplify = FALSE)
#'
#' # Do the test
#' lam <- matrix(10, 3, 3)
#' diag(lam) <- 1
#' ft <- fused.test(Ylist, Tlist, lambda = lam, n.permutations = 500)
#'
#' @export
print.ptest <- function(x, digits = 4L, ...) {
  x$n.extreme <- sum(x$null.dist >= x$observed)
  x$n.permutations <- length(x$null.dist)
  x$p.val.unbiased <- x$n.extreme/x$n.permutations
  x$p.val.biased <- (x$n.extreme + 1)/(x$n.permutations + 1)
  pval <- format.pval(x$p.val.unbiased, digits = digits,
                      eps = 1/x$n.permutations)

  cat("\nScore-based permutation test\n\n")
  cat("Null hypothesis: Population precision matrices are equal\n")
  cat("Alternative:     Population precision matrices are not equal\n\n")
  cat(sprintf("Observed statistic: U = %0.3f, ", x$observed))
  cat(sprintf("p-value %s\n", ifelse(grepl("<", pval), pval, paste("=", pval))))
  cat("Summary of null distribution obtained by permutation:\n")
  print(summary(x$null.dist, digits = digits, ...))
  return(invisible(x))
}


#' @rdname print.ptest
#' @export
summary.ptest <- function(object, ...) {
  object <- print.ptest(object, ...)
  cat("\nThe number of extreme observations under the null hypothesis")
  cat(sprintf("\nwas %d out of %d permutations.",
              object$n.extreme, object$n.permutations))

  return(invisible(object))
}


#' @rdname plot.ptest
#' @export
hist.ptest <- function(x, add.extra = TRUE, ...) {
  hist.args <- list(...)
  if (!hasArg("xlim")) {
    hist.args$xlim <- range(x$null.dist, x$observed)
  }
  if (!hasArg("col")) {
    hist.args$col <- "gray"
  }
  if (!hasArg("main")) {
    hist.args$main <- "Null distribution of U"
  }
  if (!hasArg("xlab")) {
    hist.args$xlab <- "U"
  }
  out <- do.call(hist, c(list(x$null.dist), hist.args), quote = FALSE)
  out$xname <- "x$null.dist"

  if (add.extra) {
    rug(x$null.dist)
    abline(v = x$observed, col = "red", lwd = 2)
    text(x$observed, y = par()$usr[4],
         labels = "Observed U", pos = 3, xpd = TRUE)
  }
  return(invisible(out))
}



#' Plot the results of a fusion test
#'
#' Plot a histogram of the null distribution and the observed test statistic in
#' a permutation type "fusion test".
#'
#' \code{plot.ptest} is simply a wrapper for \code{hist.ptest}.
#'
#' @param x A \code{ptest} object (a list). Usually the output of
#'   \code{\link{fused.test}}.
#' @param add.extra A logical. Add extra information to the plot.
#' @param \dots Arguments passed to \code{plot}.
#'
#' @return Invisibly returns \code{x} with extra additions.
#'
#' @author Anders Ellern Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>,
#'   Wessel N. van Wieringen
#'
#' @seealso \code{\link{fused.test}}, \code{\link{print.ptest}}
#'
#' @references Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., and
#'   van Wieringen, W.N. (2020).  Targeted Fused Ridge Estimation of Inverse
#'   Covariance Matrices from Multiple High-Dimensional Data Classes.  Journal
#'   of Machine Learning Research, 21(26): 1-52.
#'
#' @examples
#' ns <- c(10, 5, 23)
#' Ylist <- createS(ns, p = 15, topology = "banded", dataset = TRUE)
#'
#' # Use the identity target matrix for each class
#' Tlist <- replicate(length(ns), diag(15), simplify = FALSE)
#'
#' # Do the test
#' lam <- matrix(10, 3, 3)
#' diag(lam) <- 1
#' ft <- fused.test(Ylist, Tlist, lambda = lam, n.permutations = 500)
#'
#' # The returned object can alo be plotted via
#' hist(ft)
#' # or via the alias
#' plot(ft)
#' @export
plot.ptest <- function(x, add.extra = TRUE, ...) {
  hist.ptest(x, add.extra = add.extra, ...)
}



##------------------------------------------------------------------------------
##
## Sparsification and network stats
##
##------------------------------------------------------------------------------


#' Determine support of multiple partial correlation/precision matrices
#'
#' A simple wrapper for \code{\link{sparsify}} which determines the support of
#' a \code{list} of partial correlation/precision matrix by various methods and
#' returns the sparsified matrices.
#'
#' @param Plist A \code{list} of \code{numeric} precision matrices.
#' @param \dots Arguments passed to \code{\link{sparsify}}.
#'
#' @return A \code{list} of the same length as \code{Plist} with the output from
#'   \code{\link{sparsify}}.
#'
#' @author Anders Ellern Bilgrau, Wessel N. van Wierigen, Carel F.W. Peeters
#'   <cf.peeters@@vumc.nl>
#'
#' @seealso \code{\link{sparsify}}
#'
#' @examples
#' ns <- c(10, 11)
#' Ylist <- createS(ns, p = 16, dataset = TRUE)
#' Slist <- lapply(Ylist, covML)
#' Tlist <- default.target.fused(Slist, ns)
#'
#' # Obtain regularized precision under optimal penalty
#' opt <- optPenalty.fused.auto(Ylist, Tlist, cv.method = "aLOOCV",
#'                             maxit.ridgeP.fused = 1500)
#' # Use the optimal penalties
#' Plist <- ridgeP.fused(Slist, ns, lambda = opt$lambda, maxit = 1000)
#'
#' # Determine support regularized (standardized) precision under optimal penalty
#' res <- sparsify.fused(Plist, threshold = "top", verbose = FALSE)
#' round(res[[1]]$sparsePrecision, 1)
#' round(res[[2]]$sparsePrecision, 1)
#'
#' @export sparsify.fused
sparsify.fused <- function(Plist, ...) {
  return(lapply(Plist, sparsify, ...))
}



#' Gaussian graphical model network statistics
#'
#' Compute various network statistics from a \code{list} sparse precision
#' matrices. The sparse precision matrix is taken to represent the conditional
#' independence graph of a Gaussian graphical model. This function is a simple
#' wrapper for \code{\link{GGMnetworkStats}}.
#'
#' For details on the columns see \code{\link{GGMnetworkStats}}.
#'
#' @param Plist A \code{list} of sparse precision/partial correlation matrix.
#'
#' @return A \code{data.frame} of the various network statistics for each
#' class. The names of \code{Plist} is prefixed to column-names.
#'
#' @author Anders E. Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel
#' N. van Wieringen
#'
#' @seealso \code{\link{GGMnetworkStats}}
#'
#' @examples
#' ## Create some "high-dimensional" data
#' set.seed(1)
#' p <- 10
#' ns <- c(5, 6)
#' Slist <- createS(ns, p)
#'
#' ## Obtain sparsified partial correlation matrix
#' Plist    <- ridgeP.fused(Slist, ns, lambda = c(5.2, 1.3), verbose = FALSE)
#' PCsparse <- sparsify.fused(Plist , threshold = "absValue", absValueCut = 0.2)
#' SPlist <- lapply(PCsparse, "[[", "sparsePrecision") # Get sparse precisions
#'
#' ## Calculate GGM network statistics in each class
#' \dontrun{GGMnetworkStats.fused(SPlist)}
#'
#' @export GGMnetworkStats.fused
GGMnetworkStats.fused <- function(Plist) {
  res <- lapply(Plist, GGMnetworkStats, as.table = TRUE)
  if (is.null(names(res))) {
    names(res) <- seq_along(Plist)
  }
  return(as.data.frame(res))
}



#' Fused gaussian graphical model node pair path statistics
#'
#' A simple wrapper for \code{\link{GGMpathStats}}.
#'
#' @param sparsePlist A \code{list} of sparsified precision matrices.
#' @param \dots Arguments passed to \code{\link{GGMpathStats}}.
#'
#' @return A \code{list} of path stats.
#'
#' @note The function currently fails if no paths are present in one of the
#'   groups.
#'
#' @author Anders E. Bilgrau, Carel F.W. Peeters <cf.peeters@@vumc.nl>, Wessel
#'   N. van Wieringen
#'
#' @seealso \code{\link{GGMpathStats}}
#'
#' @examples
#' ## Obtain some (high-dimensional) data
#' set.seed(1)
#' ns <- c(10, 11)
#' Slist <- createS(ns, p = 7, topology = "banded")
#' Tlist <- default.target.fused(Slist, ns)
#'
#' ## Obtain regularized precision and sparsify
#' Plist <- ridgeP.fused(Slist, ns, Tlist, lambda = c(1, 1.6))
#' sparsePlist <- sparsify.fused(Plist, threshold = "absValue", absValueCut = 0.20)
#' SPlist <- lapply(sparsePlist, "[[", "sparsePrecision")
#'
#' ## Obtain information on mediating and moderating paths between nodes 14 and 23
#' res <- GGMpathStats.fused(SPlist, node1 = 3, node2 = 4, graph = FALSE)
#'
#' @export GGMpathStats.fused
GGMpathStats.fused <- function(sparsePlist, ...) {
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

