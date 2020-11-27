.cvlPchordal <- function(lambda,
                         Y,
			 target=default.target(covML(Y)),
			 zeros, cliques,
			 separators,
			 type="Alt"){
	slh <- numeric()
	for (i in 1:nrow(Y)){
		S   <- covML(Y[-i, ])
		slh <- c(slh, .LL(crossprod(Y[i, , drop=FALSE]),
		                  ridgePchordal(S,
				                lambda,
						target=target,
						zeros=zeros,
						cliques=cliques,
						separators=separators,
						type=type,
						verbose=FALSE)))
	}
	return(mean(slh))
}



.zeros2chordalCIG <- function(zeros, nCovariates){
	####################################################################
	# converts CI pattern of error precision (specified as zeros)
	# into graph-object, which is triangulated and decomposed
	####################################################################

	# zeros to adjacency matrix
	adjMat <- matrix(1, nCovariates, nCovariates)
	adjMat[zeros] <- 0
	adjMat[cbind(zeros[,2], zeros[,1])] <- 0
	diag(adjMat) <- 0

    	# convert adjacency into graphNel object
	G <- graph.adjacency(adjMat, mode="undirected")
	G <- igraph.to.graphNEL(G)

	# is graph chordal?
	if (!is.triangulated(G)){ G <- triangulate(G) }

        # decompose in cliques and separators
        decomp <- rip(G)

	# return decomposition
	return(decomp)
}









#' Support of the adjacency matrix to cliques and separators.
#'
#' Convert the support of an undirected, chordal graph into a lists of cliques
#' and separators. When the graph is not chordal, it is triangulated to make it
#' so. The undirected graph may be specified as an adjacency matrix, or by the
#' complement of its support as a matrix with the indices of the adjancency
#' matrix corresponding to absent edges. The function thus caters for the two
#' different types of output from the
#' \code{\link[rags2ridges:sparsify]{sparsify}}-function. The function is meant
#' to preceede the \code{\link{ridgePchordal}}, as it its output directly feeds
#' into the latter.
#'
#' Essentially, it is a wrapper for the \code{rip}-function from the
#' \code{gRbase}-package, which takes different input and yields slightly
#' different output. Its main purpose is to mold the input such that it is
#' convenient for the \code{ridgePchordal}-function, which provides ridge
#' maximum likelihood estimation of the precision matrix with known support.
#'
#' @param adjMat Adjacency matrix of an undirected graph.
#' @param nNodes Positive \code{integer} of length one: number nodes of the
#' network.
#' @param zeros A \code{matrix} with indices of entries of the adjacency matrix
#' that are zero. The matrix comprises two columns, each row corresponding to
#' an entry of the adjacency matrix.
#' @param verbose A \code{logical} indicator: should intermediate output be
#' printed on the screen?
#' @return A \code{list}-object comprising three slots: 'zeros', 'cliques,
#' 'separators' and 'addedEdges'. The 'zeros'-slot: a \code{matrix} with
#' indices of entries of the adjacency matrix that are zero. The matrix
#' comprises two columns, each row corresponding to an entry of the adjacency
#' matrix. The first column contains the row indices and the second the column
#' indices. The specified graph should be undirected and decomposable. If not,
#' it is symmetrized and triangulated. Hence, it may differ from the input
#' 'zeros'. The 'cliques'-slot: a \code{list}-object containing the node
#' indices per clique as obtained from the \code{rip}-function. The
#' 'separators'-slot: a \code{list}-object containing the node indices per
#' clique as obtained from the \code{rip}-function. The 'addedEdges'-slot: a
#' \code{matrix} with indices of edges that have been added in the
#' triangulation.
#' @author Wessel N. van Wieringen.
#' @seealso \code{\link[rags2ridges:sparsify]{sparsify}},
#' \code{\link{ridgePchordal}}, \code{\link[gRbase:graph-rip]{gRbase::rip}}.
#' @references Lauritzen, S.L. (2004). \emph{Graphical Models}. Oxford
#' University Press.
#' @examples
#'
#' # obtain some (high-dimensional) data
#' p <- 8
#' n <- 100
#' set.seed(333)
#' Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'
#' # create sparse precision
#' P <- covML(Y)
#' P[1:3, 6:8] <- 0
#' P[6:8, 1:3] <- 0
#'
#' # draw some data
#' S <- covML(matrix(rnorm(n*p), nrow = n, ncol = p))
#'
#' # obtain (triangulated) support info
#' zeros <- which(P==0, arr.ind=TRUE)
#' supportP <- support4ridgeP(adjMat=adjacentMat(P))
#'
#' # alternative specification of the support
#' zeros <- which(P==0, arr.ind=TRUE)
#' supportP <- support4ridgeP(nNodes=p, zeros=zeros)
#'
#' # estimate precision matrix with known (triangulated) support
#' Phat <- ridgePchordal(S, 0.1, zeros=supportP$zeros,
#' 	cliques=supportP$cliques, separators=supportP$separators)
#'
#' @export
support4ridgeP <- function(adjMat=NULL,
                           nNodes=NULL,
                           zeros=NULL,
                           verbose=FALSE){

	########################################################################
	#
	# DESCRIPTION:
	# The zero pattern of an adjancency matrix is converted into its
	# junction tree (as a list of its cliques and separators). When a
	# graph with nondecomposable support is provided, it is triangulated.
	#
	# ARGUMENTS:
	# -> adjMat     : Adjacency matrix of an undirected graph.
	# -> nNodes     : Number of nodes of the graph.
	# -> zeros      : Matrix with indices of entries of the adjacency
	#                 matrix that are constrained to zero. The matrix
	#                 comprises two columns, each row corresponding to
	#                 an entry of the adjacency matrix. The first
	#                 column contains the row indices and the second
	#                 the column indices. The specified graph should be
	#                 undirected and decomposable. If not, it is
	#                 symmetrized and triangulated.
	# -> verbose	: Logical indicating whether progress should be reported.
	#
	#
	# DEPENDENCIES:
	# require("igraph")          # functions from package : graph.adjancency,
	#                                                       igraph.to.graphNEL
	# require("gRbase")          # functions from package : triangulate
	# require("RGBL")            # functions from package : is.triangulated
	# require("graph")           # functions from package : 'graphAM'-class
	#
	########################################################################

	# iput checks
	if (is.null(adjMat) & is.null(nNodes) & is.null(zeros)){
		stop("Support not sufficiently specified.")
	}
	if (is.null(adjMat) & is.null(nNodes)){
		stop("Support not sufficiently specified.")
	}
	if (is.null(adjMat) & is.null(zeros)){
		stop("Support not sufficiently specified.")
	}
	if (!is.null(adjMat)){
		if (!is.matrix(adjMat)){
			stop("adjMat should be a matrix.")
		}
		if (nrow(adjMat) != ncol(adjMat)){
			stop("adjMat should be square matrix.")
		}
		if (!all(sort(unique(as.numeric(adjMat))) == c(0,1))){
			stop("Elements of adjMat ill-specified.")
		}
		if (!isSymmetric(adjMat)){
			stop("adjMat does not correspond to an undirect graph.")
		}
		if (!all(diag(adjMat) == 0)){
			stop("adjMat contains nonzero diagonal element.")
		}
	}
	if (!is.null(nNodes)){
		if (class(nNodes) != "numeric" & class(nNodes) != "integer"){
			stop("Input (nNodes) is of wrong class")
		}
		if (length(nNodes) != 1){
			stop("Input (nNodes) must be a scalar")
		}
	      	if (!.is.int(nNodes)){
			stop("Input (nNodes) should be an integer")
		}
		if (nNodes <= 0){
			stop("Input (nNodes) must be strictly positive")
		}
	}
	if (!is.null(zeros)){
		if (!is.null(zeros) & !inherits(zeros, "matrix")){
			stop("Input (zeros) is of wrong class.")
		}
		if (!is.null(zeros)){
			if(ncol(zeros) != 2){
				stop("Wrong dimensions of the (zeros).")
			}
		}
		if (!is.null(zeros)){
			zeros <- zeros[order(zeros[,2], zeros[,1]),]
		}
	}
	if (as.character(class(verbose)) != "logical"){
		stop("Input (verbose) is of wrong class.")
	}

	# convert zero pattern into adjacency matrix
	if (!is.null(zeros) & !is.null(nNodes) & is.null(adjMat)){
		adjMat                              <- matrix(1, nNodes, nNodes)
		adjMat[zeros]                       <- 0
		adjMat[cbind(zeros[,2], zeros[,1])] <- 0
		diag(adjMat)                        <- 0
	}

	# convert adjacency into graphNel object
	G <- igraph::igraph.to.graphNEL(igraph::graph.adjacency(adjMat,
	                                                        mode="undirected"))

	# is graph complete?
	is.complete(G)

	# check whether a chordal support has been provided
	addedEdges <- matrix(nrow=0, ncol=2)
	if (!RBGL::is.triangulated(G)){
		if(verbose){
			cat("-> provided zero pattern not chordal   : support is triangulated,", "\n")
		}
		nEdgesOld    <- numEdges(G)
		G            <- gRbase::triangulate(G)
		nEdgesNew    <- numEdges(G)
		addedEdges   <- which(as(G, "graphAM")@adjMat - adjMat == 1,
		                      arr.ind=TRUE)
		adjMat       <- as(G, "graphAM")@adjMat
		diag(adjMat) <- 1
		zeros        <- which(adjMat == 0, arr.ind=TRUE)
		if(verbose){
			cat(paste("->                                       ",
			          nEdgesNew - nEdgesOld,
			          " extra nonzeros", sep=""), "\n")
		}
	}
	if (is.null(zeros)){
		diag(adjMat) <- 1
		zeros        <- which(adjMat == 0, arr.ind=TRUE)
	}

	# decompose in cliques and separators
	decomp     <- gRbase::rip(G)
	Cliques    <- decomp$cliques
	Separators <- decomp$separators

	# convert index vectors from character to integer
	for (k in 1:length(Cliques)){
		Cliques[[k]] <- as.integer(Cliques[[k]])
	}
	for (k in 1:length(Separators)){
		Separators[[k]] <- as.integer(Separators[[k]])
	}

	return(list(zeros=zeros,
	            cliques=Cliques,
	            separators=Separators,
	            addedEdges=addedEdges))
}








#' Ridge estimation for high-dimensional precision matrices with known chordal
#' support
#'
#' Calculates various ridge estimators for high-dimensional
#' precision matrices with known support. This support should form a chordal
#' graph. If the provided support is not chordal, the function makes it so.
#'
#' @details
#' Sister function to the \code{\link{ridgeP}}-function, incorporating a chordal
#' zero structure of the precision matrix.
#'
#' The loss function for \code{type="ArchII"} is:
#'   \deqn{
#'     \log(| \mathbf{\Omega} |) - \mbox{tr}( \mathbf{S} \mathbf{\Omega} ) +
#'     \lambda
#'     \big\{
#'       \log(| \mathbf{\Omega} |) -
#'       \mbox{tr}[ (\mathbf{S} + (1+\lambda) \mathbf{I}_{p \times p})
#'       \mathbf{\Omega} ]
#'     \big\}.
#'   }
#' For \code{type="ArchI"} it is:
#'   \deqn{
#'     (1-\lambda)
#'     \big[ \log(| \mathbf{\Omega} |) -
#'           \mbox{tr}( \mathbf{S} \mathbf{\Omega} ) \big] +
#'     \lambda \big[ \log(| \mathbf{\Omega} |) - \mbox{tr}( \mathbf{\Omega} )
#'     \big],
#'   }
#' which is obtained from:
#'   \deqn{
#'     \log(| \mathbf{\Omega} |) - \mbox{tr}( \mathbf{S} \mathbf{\Omega} ) +
#'     \nu \big[ \log(| \mathbf{\Omega} |) - \mbox{tr}( \mathbf{\Omega} ) \big]
#'   }
#' by division of \eqn{(1+\nu)} and writing \eqn{\lambda = \nu / (1 + \nu)}.
#'
#' An explicit expression for the minimizer of the loss functions implied by the
#' archetypal ridge estimators (\code{type="ArchI"} and \code{type="ArchII"})
#' exists. For the simple case in which the graph decomposes into cliques
#' \eqn{\mathcal{C}_1}, \eqn{\mathcal{C}_2} and separator \eqn{\mathcal{S}} the
#' estimator is:
#'   \deqn{
#'     \widehat{\mathbf{\Omega}}  =
#'     \left(
#'       \begin{array}{lll}
#'       \, [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{C}_1 \setminus \mathcal{S}, \mathcal{C}_1 \setminus \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{C}_1 \setminus \mathcal{S}, \mathcal{S}} & \mathbf{0}_{|\mathcal{C}_1 \setminus \mathcal{S}| \times |\mathcal{C}_2 \setminus \mathcal{S}|}
#'       \\
#'       \, [\widehat{\mathbf{\Omega}}^{(\mathcal{C}_1)}]_{\mathcal{S}, \mathcal{C}_1 \setminus \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}]_{\mathcal{S}, \mathcal{S}} + [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{S}, \mathcal{S}} - \widehat{\mathbf{\Omega}}^{(\mathcal{S})} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{S}, \mathcal{C}_2 \setminus \mathcal{S}}
#'       \\
#'       \, \mathbf{0}_{|\mathcal{C}_2 \setminus \mathcal{S}| \times |\mathcal{C}_1 \setminus \mathcal{S}|} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{C}_2 \setminus \mathcal{S}, \mathcal{S}} & [\widehat{\mathbf{\Omega}}^{({\mathcal{C}_2})}]_{\mathcal{C}_2 \setminus \mathcal{S}, \mathcal{C}_2 \setminus \mathcal{S}}
#'       \end{array}
#'       \right),
#'   }
#' where \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}},
#' \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{C}_1})}} and
#' \eqn{\widehat{\mathbf{\Omega}}^{({\mathcal{S}})}} are the marginal ridge ML
#' covariance estimators for cliques \eqn{\mathcal{C}_1}, \eqn{\mathcal{C}_2}
#' and separator \eqn{\mathcal{S}}. The general form of the estimator,
#' implemented here, is analogous to that provided in Proposition 5.9 of
#' Lauritzen (2004). The proof that this estimator indeed optimizes the
#' corresponding loss function is fully analogous to that of Proposition 5.6 of
#' Lauritzen (2004).
#'
#' In case, \code{type="Alt"} no explicit expression of the maximizer of the
#' ridge penalized log-likelihood exists. However, an initial estimator
#' analogous to that for \code{type="ArchI"} and \code{type="ArchII"} can be
#' defined. In various boundary cases (\eqn{\lambda=0}, \eqn{\lambda=\infty},
#' and \eqn{\mathcal{S} = \emptyset}) this initial estimator actually optimizes
#' the loss function. In general, however, it does not. Nevertheless, it
#' functions as well-educated guess for any Newton-like optimization method:
#' convergence is usually achieved quickly. The Newton-like procedure optimizes
#' an unconstrained problem equivalent to that of the penalized log-likelihood
#' with known zeros for the precision matrix (see Dahl \emph{et al}., 2005 for
#' details).
#'
#' @param S Sample covariance \code{matrix}.
#' @param lambda A \code{numeric} representing the value of the penalty
#'   parameter.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#'   estimators.
#' @param zeros \code{Matrix} with indices of entries of the adjacency matrix
#'   that are zero. The matrix comprises two columns, each row corresponding to
#'   an entry of the adjacency matrix. The first column contains the row indices
#'   and the second the column indices. The specified graph should be undirected
#'   and decomposable. If not, use the \code{\link{support4ridgeP}} to
#'   symmetrize and triangulate. This is done automatically if \code{cliques}
#'   and \code{separators} arguments are empty lists (and the then employed
#'   \code{zeros}-object may differ from the one provided as input).
#' @param cliques A \code{list}-object containing the node indices per clique as
#'   obtained from the \code{\link{support4ridgeP}}-function.
#' @param separators A \code{list}-object containing the node indices per
#'   separator as obtained from the \code{\link{support4ridgeP}}-function.
#' @param type A \code{character} indicating the type of ridge estimator to be
#'   used. Must be one of: \code{Alt} (default), \code{ArchI}, \code{ArchII}.
#' @param optimizer A \code{character} (either \code{nlm} (default) or
#'   \code{optim}) specifying which optimization function should be used:
#'   \code{\link[stats:nlm]{nlm}} (default) or \code{\link[stats:optim]{optim}}?
#' @param grad A \code{logical} indicator: should, next to the precision matrix
#'   estimate, also the gradient be returned?
#' @param verbose A \code{logical} indicator: should intermediate output be
#'   printed on the screen?
#' @param ...  Additional arguments passed on to either
#'   \code{\link[stats:nlm]{nlm}} or \code{\link[stats:optim]{optim}}.
#'
#' @return
#' If \code{grad=FALSE}, the function returns a regularized precision
#' \code{matrix} with specified chordal sparsity structure.
#'
#' If \code{grad=TRUE}, a list is returned comprising of \emph{i)} the
#' estimated precision matrix, and \emph{ii)} the gradients at the initial and
#' at the optimal (if reached) value. The gradient is returned and it can be
#' checked whether it is indeed (close to) zero at the optimum.
#'
#' @author Wessel N. van Wieringen.
#'
#' @seealso \code{\link[rags2ridges:ridgeP]{ridgeP}}
#'
#' @references
#' Dahl, J., Roychowdhury, V., Vandenberghe, L. (2005), "Maximum
#' likelihood estimation of Gaussian graphical models: numerical implementation
#' and topology selection", Technical report, UCLA, 2005.
#'
#' Lauritzen, S.L. (2004). \emph{Graphical Models}. Oxford University Press.
#'
#' Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge estimation of
#' the VAR(1) model and its time series chain graph from multivariate
#' time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191.
#'
#' @examples
#' # obtain some (high-dimensional) data
#' p <- 8
#' n <- 100
#' set.seed(333)
#' Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'
#' # define zero structure
#' S <- covML(Y)
#' S[1:3, 6:8] <- 0
#' S[6:8, 1:3] <- 0
#' zeros <- which(S==0, arr.ind=TRUE)
#'
#' # obtain (triangulated) support info
#' supportP <- support4ridgeP(nNodes=p, zeros=zeros)
#'
#' # estimate precision matrix with known (triangulated) support
#' Phat <- ridgePchordal(S, 0.1, zeros=supportP$zeros,
#' 	cliques=supportP$cliques, separators=supportP$separators)
#'
#' @importFrom stats nlm optim
#' @importFrom igraph graph.adjacency igraph.to.graphNEL
#' @importFrom gRbase triangulate
#' @importFrom RBGL is.triangulated
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
ridgePchordal <- function(S,
                          lambda,
                          zeros,
                          cliques = list(),
                          separators = list(),
                          target = default.target(S),
                          type = "Alt",
                          optimizer = "nlm",
                          grad = FALSE,
                          verbose = TRUE,
                          ...){

	# TODO: Currently, uses the full inverse instead of the partial inverse.
	#    To be fixed in the far future.

	# Input checks
	if (!is.matrix(S)){
		stop("S should be a matrix")
	}
	if (!isSymmetric(S)){
		stop("S should be a symmetric matrix")
	}
	if (length(lambda) != 1){
		stop("Input (lambda) is of wrong length.")
	}
	if (is.na(lambda)){
		stop("Input (lambda) is not a positive number.")
	}
	if (lambda <= 0){
		stop("Input (lambda) is not a positive number.")
	}
	if (!is.null(zeros) & !inherits(zeros, "matrix")){
		stop("Input (zeros) is of wrong class.")
	}
	if (!is.null(zeros)){
		if(ncol(zeros) != 2){
			stop("Wrong dimensions of the (zeros).")
		}
	}
	if (!is.null(zeros)){
		zeros <- zeros[order(zeros[,2], zeros[,1]),]
	}
	if (!(type %in% c("Alt", "ArchI", "ArchII"))) {
		stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
	}
	if (!(optimizer %in% c("nlm", "optim"))){
		stop("type should be one of {'nlm', 'optim'}")
	}
	if (as.character(class(verbose)) != "logical"){
		stop("Input (verbose) is of wrong class.")
	}
	if (!is.list(cliques)){
		stop("Input (cliques) is of wrong class.")
	}
	if (!is.list(separators)){
		stop("Input (separators) is of wrong class.")
	}

	# intermediate output
	if(verbose){
		cat(" ", "\n")
		cat("Progress report ....", "\n")
		cat(paste("-> ----------------------------------------------------------------------",
		          sep=""), "\n")
	}

	# if chordal decomposition not supplied as a clique and separator list, make it so
	if (length(cliques) == 0){
		supportInfo <- support4ridgeP(nNodes=nrow(S), zeros=zeros);
		cliques     <- supportInfo$cliques;
		separators  <- supportInfo$separators;
	}

	# make adjacency matrix of support
	Pinit <- matrix(1, ncol=ncol(S), nrow=nrow(S));
	Pinit[zeros] <- 0;
	diag(Pinit) <- 1;

	# ensure the target has same support as chordal graph
	target[zeros] <- 0;

	# obtain the graph components
	diag(Pinit) <- 0
	components <- igraph::clusters(igraph::graph.adjacency(Pinit, mode="undirected"))$membership

	# construct init estimate
	Pinit <- .armaRidgePchordalInit(S=S,
	                                lambda=lambda,
                                        target=target,
	                                type=type,
	                                Cliques=cliques,
	                                Separators=separators)
	evs <- eigen(Pinit, only.values=TRUE)$values

	if (type=="Alt"){
		if (any(eigen(Pinit, only.values=TRUE)$values < 0)){
			Pinit <- (Pinit + diag(diag(Pinit)))/2
		}

		# number of non-converged components
		nonconvergedComp <- 0

		# report number of to be estimated parameters
		diag(Pinit) <- diag(Pinit) / 2
		Pinit[upper.tri(Pinit)] <- 0
		nonzeros <- which(Pinit != 0, arr.ind=TRUE)
		Xinit <- Pinit[nonzeros]
		if (verbose){
			cat(paste("-> optimization over                   : ",
			          nrow(nonzeros),
			          " out of ",
			          ncol(S) * (ncol(S) + 1) / 2,
			          " unique", sep=""),
			          "\n")
			cat(paste("->                                       ",
			          "parameters (",
			          round(200 * nrow(nonzeros) /  (ncol(S) * (ncol(S) + 1)), digits=2),
			          "%)",
			          sep=""), "\n")
			cat(paste("-> cond. number of initial estimate    : ",
			          round(max(evs) / max(min(evs), 0), digits=2),
			          " (if >> 100, consider", sep=""),
			          "\n")
			cat(paste("->                                       larger values of lambda)"),
			          "\n")
			cat(paste("-> # graph components                  : ",
			          length(unique(components)),
			          " (optimization per component)", sep=""),
			          "\n")
			cat(paste("-> optimization per component          : ", sep=""),
			          "\n")
		}

		if (verbose){
			pBar <- utils::txtProgressBar(min=0,
			                              max=abs(length(unique(components))),
			                              style=3,
			                              char=".")
		}

		# estimate per graph component
		for (subG in unique(components)){
			if (verbose){
				utils::setTxtProgressBar(pBar, subG);
			  Sys.sleep(10^(-10))  # TODO: WHY?
			}

			# construct component data
			ids      <- which(components == subG)
			Psub     <- Pinit[ids, ids, drop=FALSE]
			nonzeros <- which(Psub != 0, arr.ind=TRUE)
			x0       <- Psub[nonzeros]
			rm(Psub)

			# matrices for alternative parametrization of the sparse precision matrix
			E1 <- matrix(0, ncol=nrow(nonzeros), nrow=length(ids))
			E1[cbind(nonzeros[,1], c(1:nrow(nonzeros)))] <- 1
			E2 <- matrix(0, ncol=nrow(nonzeros), nrow=length(ids))
			E2[cbind(nonzeros[,2], c(1:nrow(nonzeros)))] <- 1

			if (nrow(nonzeros) != length(ids) * (length(ids) + 1) / 2){
				# minimize minus penalized log-likelihood
				if (optimizer=="optim"){
					xhat <- optim(x0,
					              .armaPenLLreparP,
					              gr=.armaPenLLreparPgrad,
					              method="BFGS",
					              E1=E1,
					              E2=E2,
					              S=S[ids,ids,drop=FALSE],
					              lambda=lambda,
					              target=target[ids,ids,drop=FALSE],
					              nonzerosR=nonzeros[,1],
					              nonzerosC=nonzeros[,2],
					              ...)
					if(xhat$convergence != 0){
						nonconvergedComp <- nonconvergedComp + 1
					}
					Pinit[ids, ids] <- E1 %*% diag(xhat$par, ncol=length(x0)) %*% t(E2) +
					                   E2 %*% diag(xhat$par, ncol=length(x0)) %*% t(E1)
				}

				if (optimizer=="nlm"){
					xhat <- nlm(.armaPenLLreparPforNLM,
					            p=x0,
					            E1=E1,
					            E2=E2,
					            S=S[ids,ids,drop=FALSE],
					            lambda=lambda,
					            target=target[ids,ids,drop=FALSE],
					            nonzerosR=nonzeros[,1],
					            nonzerosC=nonzeros[,2],
					            check.analyticals=FALSE,
					            ...)
					if(xhat$code != 0 & xhat$code != 1){
						nonconvergedComp <- nonconvergedComp + 1
					}
					Pinit[ids, ids] <- E1 %*% diag(xhat$estimate, ncol=length(x0)) %*% t(E2) +
					                   E2 %*% diag(xhat$estimate, ncol=length(x0)) %*% t(E1)
				}
			} else {
				Pinit[ids, ids] <- E1 %*% diag(x0, ncol=length(x0)) %*% t(E2) +
				                   E2 %*% diag(x0, ncol=length(x0)) %*% t(E1)
			}
		}
	}

	if (verbose){
		if (type=="Alt"){
			cat("\n");
			cat("\n")
		}
		cat(paste("-> estimation done ...", sep=""), "\n")
		cat(paste("-> formatting output ...", sep=""), "\n")
	}

	# reformatting for reporting the gradient
	diag(Pinit)             <- diag(Pinit) / 2
	Pinit[upper.tri(Pinit)] <- 0
	nonzeros                <- which(Pinit != 0, arr.ind=TRUE)
	x0                      <- Pinit[nonzeros]
	E1                      <- matrix(0, ncol=nrow(nonzeros), nrow=nrow(Pinit))
	E1[cbind(nonzeros[,1], c(1:nrow(nonzeros)))] <- 1
	E2                      <- matrix(0, ncol=nrow(nonzeros), nrow=nrow(Pinit))
	E2[cbind(nonzeros[,2], c(1:nrow(nonzeros)))] <- 1
	rm(Pinit)

	if (verbose & type=="Alt"){
		cat(paste("-> overall summary ....", sep=""), "\n")
		cat(paste("-> initial pen. log-likelihood         : ",
			  round(- .armaPenLLreparP(Xinit, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]), 8), sep=""),
		          "\n")
		cat(paste("-> optimized pen. log-likelihood       : ",
		          round(- .armaPenLLreparP(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2]), 8), sep=""),
		          "\n")
		if (nonconvergedComp == 0){
			cat("-> optimization                        : converged (most likely)", "\n")
			cat("->                                       for all components", "\n")
		}
		if (nonconvergedComp > 0){
			cat("-> optimization                        : for ", nonconvergedComp, " components", "\n")
			cat("->                                       max. no. iterations reached (or", "\n")
			cat("->                                       other indication of possible", "\n")
			cat("->                                       convergence failure)", "\n")
		}
	}
	if (verbose){ cat(paste("-> ----------------------------------------------------------------------", sep=""), "\n")  }

	# return precision matrix from alternative parametrization
	if (!grad){
		return(E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1))
	} else {
		if (type=="Alt"){
			return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1),
			            grad=.armaPenLLreparPgrad(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2])))
		}
		if (type=="ArchI"){
			return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1),
			            grad=.armaPenLLreparGradArchI(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2])))
		}
		if (type=="ArchII"){
			return(list(P = E1 %*% diag(x0) %*% t(E2) + E2 %*% diag(x0) %*% t(E1),
			            grad=.armaPenLLreparGradArchII(x0, E1, E2, S, lambda, target, nonzeros[,1], nonzeros[,2])))
		}
	}
}










#' Automatic search for penalty parameter of ridge precision estimator with
#' known chordal support
#'
#' Automic search for the optimal ridge penalty parameter for the ridge
#' estimator of the precision matrix with known chordal support. Optimal in the
#' sense that it yields the maximum cross-validated likelihood. The search
#' employs the Brent algorithm as implemented in the
#' \code{\link[stats:optim]{optim}} function.
#'
#' See the function \code{\link[stats:optim]{optim}} for details on the
#' implementation of the Brent algorithm.
#'
#' @param Y Data \code{matrix}. Variables assumed to be represented by columns.
#' @param lambdaMin A \code{numeric} giving the minimum value for the penalty
#' parameter.
#' @param lambdaMax A \code{numeric} giving the maximum value for the penalty
#' parameter.
#' @param lambdaInit A \code{numeric} giving the initial value for the penalty
#' parameter.
#' @param target A target \code{matrix} (in precision terms) for Type I ridge
#' estimators.
#' @param zeros A \code{matrix} with indices of entries of the precision matrix
#' that are constrained to zero. The matrix comprises two columns, each row
#' corresponding to an entry of the precision matrix. The first column contains
#' the row indices and the second the column indices. The specified conditional
#' independence graph implied by the zero-structure of the precision should be
#' undirected and decomposable. If not, it is symmetrized and triangulated.
#' @param cliques A \code{list}-object containing the node indices per clique
#' as obtained from the \code{\link{support4ridgeP}}-function.
#' @param separators A \code{list}-object containing the node indices per
#' separator as obtained from the \code{\link{support4ridgeP}}-function.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: \code{Alt}, \code{ArchI}, \code{ArchII}.
#' @return A \code{numeric} with the LOOCV optimal choice for the ridge penalty
#' parameter.
#' @author Wessel N. van Wieringen.
#' @seealso \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}},
#' \code{\link[rags2ridges:ridgeP]{ridgeP}},
#' \code{\link[rags2ridges:optPenalty.aLOOCV]{optPenalty.aLOOCV}},
#' \code{\link[rags2ridges:optPenalty.kCV]{optPenalty.kCV}}
#' @references Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge
#' estimation of the VAR(1) model and its time series chain graph from
#' multivariate time-course omics data", \emph{Biometrical Journal}, 59(1),
#' 172-191.
#'
#' Van Wieringen, W.N. and Peeters, C.F.W. (2016), "Ridge Estimation of Inverse
#' Covariance Matrices from High-Dimensional Data", \emph{Computational
#' Statistics and Data Analysis}, 103, 284-303.
#' @examples
#'
#' # generate data
#' p <- 8
#' n <- 100
#' set.seed(333)
#' Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'
#' # define zero structure
#' S <- covML(Y)
#' S[1:3, 6:8] <- 0
#' S[6:8, 1:3] <- 0
#' zeros <- which(S==0, arr.ind=TRUE)
#'
#' # obtain (triangulated) support info
#' supportP <- support4ridgeP(nNodes=p, zeros=zeros)
#'
#' # determine optimal penalty parameter
#' \dontrun{
#' optLambda <- optPenaltyPchordal(Y, 10^(-10), 10, 0.1, zeros=supportP$zeros,
#' 	cliques=supportP$cliques, separators=supportP$separators)
#' }
#' optLambda <- 0.1
#'
#' # estimate precision matrix with known (triangulated) support
#' Phat <- ridgePchordal(S, optLambda, zeros=supportP$zeros,
#' 	cliques=supportP$cliques, separators=supportP$separators)
#'
#' @importFrom stats optim
#' @export
optPenaltyPchordal <- function (Y,
                                lambdaMin,
                                lambdaMax,
                                lambdaInit=(lambdaMin+lambdaMax)/2,
                                zeros,
                                cliques=list(),
                                separators=list(),
                                target=default.target(covML(Y)),
                                type="Alt"){

	########################################################################
	# determines the optimal value of the penalty parameter by application
	# of the Brent algorithm to the (leave-one-out) cross-validated
	# log-likelihood
	########################################################################

	# input checks
	if (!inherits(Y, "matrix")){
		stop("Input (Y) is of wrong class.")
	}
	if (sum(is.na(Y)) != 0){
		stop("Matrix Y contains missings.")
	}
	if (as.character(class(lambdaMin)) != "numeric"){
		stop("Input (lambdaMin) is of wrong class")
	}
	if (length(lambdaMin) != 1){
		stop("lambdaMin must be a scalar")
	}
	if (lambdaMin <= 0){
		stop("lambdaMin must be positive")
	}
	if (class(lambdaMax) != "numeric"){
		stop("Input (lambdaMax) is of wrong class")
	}
	if (length(lambdaMax) != 1){
		stop("lambdaMax must be a scalar")
	}
	if (lambdaMax <= lambdaMin){
		stop("lambdaMax must be larger than lambdaMin")
	}
	if (as.character(class(lambdaInit)) != "numeric"){
		stop("Input (lambdaInit) is of wrong class")
	}
	if (length(lambdaInit) != 1){
		stop("lambdaInit must be a scalar")
	}
	if (lambdaInit <= lambdaMin){
		stop("lambdaInit must be larger than lambdaMin")
	}
	if (lambdaMax <= lambdaInit){
		stop("lambdaInit must be smaller than lambdaMax")
	}
	if (!(type %in% c("Alt", "ArchI", "ArchII"))){
		stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
	}
	if (!is.list(cliques)){
		stop("Input (cliques) is of wrong class.")
	}
	if (!is.list(separators)){
		stop("Input (separators) is of wrong class.")
	}

	# if chordal decomposition not supplied as a clique and separator list, make it so
	if (length(cliques) == 0){
		supportInfo <- support4ridgeP(nNodes=ncol(Y), zeros=zeros);
		cliques     <- supportInfo$cliques;
		separators  <- supportInfo$separators;
	}

	# determine optimal value of ridge penalty parameter
	optLambda <- optim(lambdaInit,
	                   .cvlPchordal,
	                   method="Brent",
	                   lower=lambdaMin,
	                   upper=lambdaMax,
	                   Y=Y,
	                   target=target,
	                   zeros=zeros,
	                   cliques=cliques,
	                   separators=separators,
	                   type=type)$par
	return(optLambda)
}










#' Ridge estimation for high-dimensional precision matrices with known sign of
#' off-diagonal precision elements.
#'
#' Function that calculates the ridge estimators for high-dimensional precision
#' matrices with known sign of the off-diagonal precision elements.
#'
#' Modified version of the \code{\link{ridgePchordal}}-function, now the ridge
#' precision matrix estimate has off-diagonal elements equalling zero or of the
#' specified sign. The estimate is found by solving a constrained estimation
#' problem. This is done numerically and employs the
#' \code{\link[stats:nlminb]{nlminb}} and
#' \code{\link[stats:constrOptim]{constrOptim}} procedure of R. These
#' procedures are initiated by the ridge ML precision estimate and its
#' off-diagonal elements with the excluded sign set to (effectively) zero.
#'
#' @param S Sample covariance \code{matrix}.
#' @param lambda A \code{numeric} representing the value of the penalty
#' parameter.
#' @param sign A character indicating the required sign of the off-diagonal
#' elements of ridge precision estimate. Must be either: "pos" (positive) and
#' "neg" (negative).
#' @param target A target \code{matrix} (in precision terms) for the ridge
#' precision estimator.
#' @param type A \code{character} indicating the type of ridge estimator to be
#' used. Must be one of: "Alt", "ArchI", "ArchII".
#' @param method A \code{character} : which optimization function should be
#' used: \code{"nlm"} (default) or \code{"optim"} which refer to
#' \code{\link[stats:nlminb]{nlminb}} or
#' \code{\link[stats:constrOptim]{constrOptim}}, respectively.
#' @param verbose \code{Logical} indicator: should intermediate output be
#' printed on the screen?
#' @param ...  Additional arguments passed on to either
#' \code{\link[stats:nlminb]{nlminb}} or
#' \code{\link[stats:constrOptim]{constrOptim}}.
#' @return The function returns a regularized precision \code{matrix} with
#' off-diagonal elements of specified signed or zero.
#' @author W.N. van Wieringen.
#' @seealso \code{\link[rags2ridges:ridgeP]{ridgeP}},
#' \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}
#' @examples
#'
#' # obtain some data
#' p <- 8
#' n <- 100
#' set.seed(333)
#' Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'
#' # obtain regularized precision matrix with off-diagonal elements of specified signed
#' ridgePsign(covML(Y), lambda=0.1, sign="pos")
#'
#' @importFrom stats constrOptim nlminb
#' @export
ridgePsign <- function(S,
                       lambda,
                       sign,
                       target = default.target(S),
                       type = "Alt",
                       method = "nlm",
                       verbose = TRUE,
                       ...){

  # Input checks
	if (!is.matrix(S)){
		stop("S should be a matrix")
	}
	if (!isSymmetric(S)){
		stop("S should be a symmetric matrix")
	}
	if (length(lambda) != 1){
		stop("Input (lambda) is of wrong length.")
	}
	if (is.na(lambda)){
		stop("Input (lambda) is not a positive number.")
	}
	if (lambda <= 0){
		stop("Input (lambda) is not a positive number.")
	}
	if (!(sign %in% c("neg", "pos"))) {
		stop("type should be one of {'neg', 'pos'}")
	}
	if (!(type %in% c("Alt", "ArchI", "ArchII"))){
		stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
	}
	if (!(method %in% c("nlm", "optim"))){
		stop("type should be one of {'nlm', 'optim'}")
	}
	if (as.character(class(verbose)) != "logical"){
		stop("Input (verbose) is of wrong class.")
	}

	# check whether a chordal support has been provided
	Phat <- ridgeP(S, lambda, type, target)
	if (sign=="pos"){
		Phat[Phat < 0] <- 0
	}
	if (sign=="neg"){
		slh            <- diag(Phat);
		Phat[Phat > 0] <- 0;
		diag(Phat)     <- slh
	}

	if (type=="Alt"){
		# construct component data
		Ptemp                   <- matrix(1, nrow=nrow(Phat), ncol=ncol(Phat))
            	Ptemp[upper.tri(Ptemp)] <- 0
            	nonzeros                <- which(Ptemp > 0, arr.ind=TRUE)
            	diag(Phat)              <- diag(Phat)/2
		x0                      <- Phat[nonzeros]
		rm(Ptemp)

		# matrices for alternative parametrization of the sparse precision matrix
		E1 <- matrix(0, ncol=nrow(nonzeros), nrow=ncol(Phat))
		E1[cbind(nonzeros[,1], c(1:nrow(nonzeros)))] <- 1
		E2 <- matrix(0, ncol=nrow(nonzeros), nrow=ncol(Phat))
		E2[cbind(nonzeros[,2], c(1:nrow(nonzeros)))] <- 1

		# minimize minus penalized log-likelihood
		if (method=="optim"){
		        if (sign == "pos"){
 	                       xhat <- constrOptim(theta=x0+10^(-10),
 			                           f=.armaPenLLreparP,
				                   grad=.armaPenLLreparPgrad,
				                   method="BFGS",
	                                           E1=E1,
	                                           E2=E2,
	                                           S=S,
				                   lambda=lambda,
	                                           target=target,
                                                   ui=diag(length(x0)),
                                                   ci=rep(0, length(x0)),
                                                   nonzerosR=nonzeros[,1],
                                                   nonzerosC=nonzeros[,2],
				                   ...)
	                }
	                if (sign == "neg"){
        	                ui   <- rep(-1, length(x0))
        	                ui[cumsum(c(1, ncol(Phat):2))] <- 1
        	                xhat <- constrOptim(theta=x0-10^(-10),
			                            f=.armaPenLLreparP,
                                                    grad=.armaPenLLreparPgrad,
                                                    method="BFGS",
                                                    E1=E1,
                                                    E2=E2,
						    S=S,
						    lambda=lambda,
					            target=target,
						    ui=diag(ui),
						    ci=rep(0, length(x0)),
						    nonzerosR=nonzeros[,1],
						    nonzerosC=nonzeros[,2],
						    ...)
        	        }
        	        if(xhat$convergence != 0){
				nonconvergedComp <- nonconvergedComp + 1
			}
	                Phat <- E1 %*% diag(xhat$par, ncol=length(x0)) %*% t(E2) +
			        E2 %*% diag(xhat$par, ncol=length(x0)) %*% t(E1)
		}

		if (method=="nlm"){
			if (sign == "pos"){
				xhat <- nlminb(x0,
					       .armaPenLLreparPforNLM,
					       E1=E1,
					       E2=E2,
					       S=S,
					       lambda=lambda,
					       target=target,
					       lower=rep(0, length(x0)),
					       nonzerosR=nonzeros[,1],
					       nonzerosC=nonzeros[,2],
					       ...)
			}
        	        if (sign == "neg"){
				upper <- rep(0, length(x0))
				upper[cumsum(c(1, ncol(Phat):2))] <- 10^10
    				xhat <- nlminb(x0,
					       .armaPenLLreparPforNLM,
				               E1=E1,
					       E2=E2,
					       S=S,
					       lambda=lambda,
					       target=target,
					       nonzerosR=nonzeros[,1],
					       nonzerosC=nonzeros[,2],
					       upper=upper,
					       ...)
    			}
			if(xhat$convergence != 0){
				converged <- TRUE
			}
			Phat <- E1 %*% diag(xhat$par, ncol=length(x0)) %*% t(E2) +
                                E2 %*% diag(xhat$par, ncol=length(x0)) %*% t(E1)
		}
	}
	return(Phat)
}
