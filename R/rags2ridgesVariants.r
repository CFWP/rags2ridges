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


ridgePchordal <- function(S,
                          lambda,
                          zeros,
                          cliques=list(),
                          separators=list(),
                          target=default.target(S),
                          type="Alt",
                          optimizer="nlm",
                          grad=FALSE,
                          verbose=TRUE,
                          ...){

	########################################################################
	#
	# DESCRIPTION:
	# Ridge estimation of the precision matrix with known zeros.
	# Nonzeros should form a chordal graph. If not chordal, the function
	# makes it so.
	#
	# ARGUMENTS:
	# -> S          	: Sample covariance matrix.
	# -> lambda     	: A numeric representing the value of the
	#                         penalty parameter.
	# -> target     	: A target matrix (in precision terms) for
	#                         Type I ridge estimators.
	# -> zeros	        : A matrix with indices of entries of the
	#                         precision matrix that are constrained to
	#                         zero. The matrix comprises two columns, each
	#                         row corresponding to an entry of the adjacency
	#                         matrix. The first column contains the row
	#                         indices and the second the column indices. The
	#                         specified graph should be undirected and
	#                         decomposable. If not, it is symmetrized and
	#                         triangulated. Hence, it may differ from the
	#                         input 'zeros'.
	# -> cliques            : A 'list'-object containing the node indices
	#                         per clique as obtained from the
	#                         'support4ridgeP'-function.
	# -> separators         : A 'list'-object containing the node indices
	#                         per separator as obtained from the
	#                         'support4ridgeP'-function.
	# -> type       	: A character indicating the type of ridge
	#                         estimator to be used. Must be one of:
	#                         "Alt", "ArchI", "ArchII".
	# -> optimizer     	: Which optimization function should be
	#                         used: "optim" or "nlm"?
	# -> grad       	: Logical indicator: should, next to the
	#                         precision matrix estimate, also the
	#                         gradient be returned?
	# -> verbose    	: Logical indicator: should intermediate
	#                         output be printed on the screen?
	# -> ...        	: Additional arguments passed on to either
	#                         "optim" or "nlm".
	#
	# DEPENDENCIES:
	# require("igraph")          # functions from package : graph.adjancency,
	#                                                       igraph.to.graphNEL
	# require("gRbase")          # functions from package : triangulate
	# require("RGBL")            # functions from package : is.triangulated
	# require("rags2ridges")     # functions from package : adjacentMat,
	#                                                       default.target
	# require("stats")           # functions from package : nlm, optim
	# require("utils")           # functions from package : txtProgressBar,
	#                                                       setTxtProgressBar,
	#                                                       Sys.sleep
	#
	# NOTES:
	# 1) Currently, uses the full inverse instead of the partial inverse.
	#    To be fixed in the far future.
	#
	########################################################################

   	# input checks
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
				Sys.sleep(10^(-10))
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




ridgePsign <- function(S,
                       lambda,
                       sign,
                       target=default.target(S),
                       type="Alt",
                       method="nlm",
                       verbose=TRUE,
                       ...){

	########################################################################
	#
	# DESCRIPTION:
	# Ridge estimation of the precision matrix with known sign of its
	# off-diagonal elements.
	#
	# ARGUMENTS:
	# -> S          : Sample covariance matrix.
	# -> lambda     : A numeric representing the value of the penalty
	#                 parameter.
	# -> target     : A target matrix (in precision terms) for Type I
	#                 ridge estimators.
	# -> sign       : A character indicating the required sign of the
	#                 off-diagonal elements of ridge precision estimate.
	#                 Must be either: "pos" (positive) and "neg" (negative).
	# -> type       : A character indicating the type of ridge estimator to
	#                 be used. Must be one of: "Alt", "ArchI", "ArchII".
	# -> method     : Which optimization function should be used:
	#                 "optim" (constrOptim) or "nlm" (nlminb)?
	# -> grad       : Logical indicator: should, next to the precision
	#                 matrix estimate, also the gradient be returned?
	# -> verbose    : Logical indicator: should intermediate output be
	#                 printed on the screen?
	# -> ...        : Additional arguments passed on to either "optim"
	#                 or "nlm".
	#
	# DEPENDENCIES:
	# require("rags2ridges")     # functions from package : ridgeP,
	#                                                       default.target
	# require("stats")           # functions from package : nlminb,
	#                                                       constrOptim
	#
	# NOTES:
	# ...
	#
	########################################################################

   	# input checks
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
