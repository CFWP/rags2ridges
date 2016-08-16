sparsifyGGM <- function(P, threshold = c("absValue", "connected", "localFDR", "top"),
                     absValueCut = .25, FDRcut = .9,
                     top = 10, output = "heavy", verbose = TRUE){
  ##############################################################################
  # - Function that sparsifies/determines support of a partial correlation
  #   matrix
  # - Support can be determined by absolute value thresholding or by local FDRs
  #   thresholding
  # - One can also choose to threshold based on the top X of absolute partial
  #   correlations
  # - Local FDR operates on the nonredundant non-diagonal elements of a partial
  #   correlation matrix
  # - Function is to some extent a wrapper around certain 'fdrtool' functions
  # - P           > (possibly shrunken) precision matrix
  # - threshold   > signifies type of thresholding
  # - absValueCut > cut-off for partial correlation elements selection based on
  #                 absolute value thresholding.
  #                 Only when threshold = 'absValue'. Default = .25
  # - FDRcut      > cut-off for partial correlation element selection based on
  #                 local FDR thresholding
  #                 Only when threshold = 'localFDR'. Default = .9
  # - top         > partial correlation element selection based on retainment
  #                 'top' number of absolute partial correlation.
  #                 Only when threshold = 'top'. Default = 10
  # - output      > must be one of {"heavy", "light"}, default = "heavy"
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen. Only when threshold = 'localFDR'. Default = TRUE.
  #
  # NOTES:
  # - Input (P) may be the partial correlation matrix or the standardized
  #   precision matrix. These are identical up to the signs of off-diagonal
  #   elements. Either can be used as it has no effect on the thresholding
  #   operator and the ensuing sparsified result.
  # - Input (P) may also be the unstandardized precision matrix. The function
  #   converts it to the partial correlation matrix
  # - The function evaluates if the input (P) is a partial
  #   correlation/standardized precision matrix or an unstandardized precision
  #   matrix. If the input amounts to the latter both the sparsified partial
  #   correlation matrix and the corresponding sparsified precision matrix are
  #   given as output (when output = "heavy"). Otherwise, the ouput consists of
  #   the sparsified partial correlation/standardized precision matrix.
  # - When output = "light", only the (matrix) positions of the zero and
  #   non-zero elements are returned.
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("fdrtool")

  if (!is.matrix(P)){
    stop("Input (P) should be a matrix")
  }
  else if (!isSymmetric(P)){
    stop("Input (P) should be a symmetric matrix")
  }
  else if (!evaluateS(P, verbose = FALSE)$posEigen){
    stop("Input (P) is expected to be positive definite")
  }
  else if (missing(threshold)){
    stop("Need to specify type of sparsification ('absValue' or 'localFDR' ",
         "or 'top')")
  }
  else if (!(threshold %in% c("absValue", "connected", "localFDR", "top"))){
    stop("Input (threshold) should be one of {'absValue', 'connected', 'localFDR', 'top'}")
  }
  else if (!(output %in% c("light", "heavy"))){
    stop("Input (output) should be one of {'light', 'heavy'}")
  }
  else {
    # Obtain partial correlation matrix
    if (all(length(unique(diag(P))) == 1 & unique(diag(P)) == 1)){
      stan = TRUE
      PC  <- P
    } else {
      stan = FALSE
      PC  <- pcor(P)
    }

    # Number of nonredundant elements
    NR <- (ncol(P)*(ncol(P) - 1))/2

    # Obtain sparsified matrix
    if (threshold == "top"){
      if (class(top) != "numeric"){
        stop("Input (top) is of wrong class")
      } else if (length(top) != 1){
        stop("Input (top) must be a scalar")
      } else if (!.is.int(top)){
        stop("Input (top) should be a numeric integer")
      } else if (top <= 0){
        stop("Input (top) must be strictly positive")
      } else if (top >= NR){
        stop("Input (top) must be smaller than the number of nonredundant ",
             "off-diagonal elements of the input matrix P")
      } else {
        absValueCut <- sort(abs(PC[upper.tri(PC)]),
                            decreasing = TRUE)[ceiling(top)]
        threshold   <- "absValue"
      }
    }

    # Obtain sparsified matrix
    if (threshold == "connected"){
	sumPC <- summary(abs(PC[upper.tri(PC)]))
	maxPC <- as.numeric(sumPC[6]); minPC <- as.numeric(sumPC[1])
	for (j in 1:100){
		absValueCut <- (maxPC + minPC)/2
        	PC0 <- PC
        	PC0[!(abs(PC0) >= absValueCut)] <- 0
		if (igraph::is.connected(graph.adjacency(adjacentMat(PC0), "undirected"))){ 
			minPC <- absValueCut	
		} else {
			maxPC <- absValueCut	
		}
		if (abs(absValueCut - (maxPC + minPC)/2) < 10^(-10)){ absValueCut <- minPC; break }
	}
        threshold   <- "absValue"
    }

    if (threshold == "absValue"){
      if (class(absValueCut) != "numeric"){
        stop("Input (absValueCut) is of wrong class")
      } else if (length(absValueCut) != 1){
        stop("Input (absValueCut) must be a scalar")
      } else if (absValueCut <= 0 | absValueCut >= 1){
        stop("Input (absValueCut) must be in the interval (0,1)")
      } else {
        PC0 <- PC
        PC0[!(abs(PC0) >= absValueCut)] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    if (threshold == "localFDR"){
      if (class(FDRcut) != "numeric"){
        stop("Input (FDRcut) is of wrong class")
      } else if (length(FDRcut) != 1){
        stop("Input (FDRcut) must be a scalar")
      } else if (FDRcut <= 0 | FDRcut >= 1){
        stop("Input (FDRcut) must be in the interval (0,1)")
      } else if (class(verbose) != "logical"){
        stop("Input (verbose) is of wrong class")
      } else {
        lFDRs <- 1 - fdrtool(PC[upper.tri(PC)], "correlation",
                             plot = verbose, verbose = verbose)$lfdr
        PC0   <- diag(nrow(PC))
        PC0[lower.tri(PC0)] <- 1
        zeros <- which(PC0 == 0, arr.ind = TRUE)
        zeros <- zeros[which(lFDRs <= FDRcut),]
        PC0   <- PC
        PC0[zeros] <- 0
        PC0[cbind(zeros[,2], zeros[,1])] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    # Return
    NNZ <- length(which(PC0[upper.tri(PC0)] != 0))
    cat("- Retained elements: ", NNZ, "\n")
    cat("- Corresponding to", round(NNZ/NR,4) * 100,"% of possible edges \n")
    cat(" \n")

    if (output == "heavy"){
      if (stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        return(PC0)
      }
      if (!stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        colnames(P0)  = rownames(P0)  <- colnames(P)
        return(list(sparseParCor = PC0, sparsePrecision = P0))
      }
    }
    if (output == "light"){
      return(list(zeros = which(PC0 == 0, arr.ind = TRUE),
                  nonzeros = which(PC0 != 0, arr.ind = TRUE)))
    }
  }
}

