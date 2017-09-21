################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Deprecated Collection
##
##------------------------------------------------------------------------------
################################################################################
################################################################################


conditionNumberPlot <- function(S, lambdaMin, lambdaMax, step, type = "Alt",
                                target = default.target(S), norm = "2",
                                digitLoss = FALSE, rlDist = FALSE,
                                vertical = FALSE, value, main = TRUE,
                                nOutput = FALSE, verbose = TRUE){
  ##############################################################################
  # - Function that visualizes the spectral condition number against the
  #   regularization parameter
  # - Can be used to heuristically determine the (minimal) value of the penalty
  #   parameter
  # - The ridge estimators operate by shrinking the eigenvalues
  # - This is especially the case when targets are used that lead to rotation
  #   equivariant estimators
  # - Maximum shrinkage (under rotation equivariance) implies that all
  #   eigenvalues will be equal
  # - Ratio of maximum and minimum eigenvalue of P can then function as
  #   a heuristic
  # - It's point of stabilization can give an acceptable value for the penalty
  # - The ratio boils down to the (spectral) condition number of a matrix
  # - S         > sample covariance/correlation matrix
  # - lambdaMin > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax > maximum value penalty parameter (dependent on 'type')
  # - step      > determines the coarseness in searching the grid
  #               [lambdaMin, lambdaMax].The steps on the grid are equidistant
  #               on the log scale
  # - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target    > target (precision terms) for Type I estimators,
  #               default = default.target(S)
  # - norm      > indicates the norm under which the condition number is to be
  #               estimated. Default is the L2-norm. The L1-norm can be
  #               (cheaply) approximated
  # - digitLoss > logical indicating if the approximate loss in digits of
  #               accuracy should also be plotted. Default = FALSE
  # - rlDist    > logical indicating if relative distance to set of singular
  #               matrices should also be plotted. Default = FALSE
  # - vertical  > optional argument for visualization vertical line in graph
  #               output, default = FALSE. Can be used to indicate the value of,
  #               e.g., the optimal penalty as indicated by some routine. Can be
  #               used to assess if this optimal penalty will lead to a
  #               well-conditioned estimate
  # - value     > indicates constant on which to base vertical line when
  #               vertical = TRUE
  # - main      > logical indicating if plot should contain type of estimator as
  #               main title
  # - nOutput   > logical indicating if numeric output should be given (lambdas
  #               and condition numbers)
  # - verbose   > logical indicating if intermediate output should be printed on
  #               screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")
  # require("Hmisc")
  # require("sfsmisc")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(S)){
    stop("S should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("S should be a symmetric matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("lambdaMin must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("lambdaMin must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("lambdaMax must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("lambdaMax must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("step should be integer")
  }
  else if (step <= 0){
    stop("step should be a positive integer")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else if (!isSymmetric(target)){
    stop("Shrinkage target should be symmetric")
  }
  else if (dim(target)[1] != dim(S)[1]){
    stop("S and target should be of the same dimension")
  }
  else if (type == "Alt" & !all(target == 0) &
           any(eigen(target, symmetric = TRUE, only.values = T)$values <= 0)){
    stop("When target is not a null-matrix it should be p.d.")
  }
  else if (type == "ArchI" & lambdaMax > 1){
    stop("lambda should be in (0,1] for this type of Ridge estimator")
  }
  else if (type == "ArchI" &
           any(eigen(target, symmetric = TRUE, only.values = T)$values <= 0)){
    stop("Target should be p.d.")
  }
  else if (!(norm %in% c("2", "1"))){
    stop("norm should be one of {'2', '1'}")
  }
  else if (class(digitLoss) != "logical"){
    stop("Input (digitLoss) is of wrong class")
  }
  else if (class(rlDist) != "logical"){
    stop("Input (rlDist) is of wrong class")
  }
  else if (digitLoss & rlDist){
    stop("Only one of 'digitLoss' and 'rlDist' may be TRUE")
  }
  else if (class(vertical) != "logical"){
    stop("Input (vertical) is of wrong class")
  }
  else if (class(main) != "logical"){
    stop("Input (main) is of wrong class")
  }
  else if (class(nOutput) != "logical"){
    stop("Input (nOutput) is of wrong class")
  }
  else {
    # Set preliminaries
    lambdas <- lseq(lambdaMin, lambdaMax, length = step)
    condNR  <- numeric()

    if (norm == "2"){
      # Calculate spectral condition number ridge estimate on lambda grid
      if (verbose){cat("Calculating spectral condition numbers...", "\n")}
      if (type == "Alt" & all(target == 0)){
        Spectral <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        for (k in 1:length(lambdas)){
          Eigshrink <- .armaEigShrink(Spectral, lambdas[k])
          condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
        }
      } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0)
                 & (length(unique(diag(target))) == 1)){
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        for (k in 1:length(lambdas)){
          Eigshrink <- .armaEigShrink(Spectral, lambdas[k], cons = varPhi)
          condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
        }
      } else {
        if (type == "Alt"){
          for (k in 1:length(lambdas)){
            P         <- .armaRidgePAnyTarget(S, target = target, lambdas[k])
            Eigs      <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
        if (type != "Alt"){
          for (k in 1:length(lambdas)){
            P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
            Eigs      <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
      }
    }

    if (norm == "1"){
      # Calculate approximation to condition number under 1-norm
      if (verbose){cat("Approximating condition number under 1-norm...", "\n")}
      if (type == "Alt"){
        for (k in 1:length(lambdas)){
          P         <- .armaRidgeP(S, target = target, lambdas[k])
          condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
        }
      }
      if (type != "Alt"){
        for (k in 1:length(lambdas)){
          P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
          condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
        }
      }
    }

    # Visualization
    if (verbose){cat("Plotting...", "\n")}
    if (main){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
    }
    if (!main){Main = " "}
    if (norm == "2"){Ylab = "spectral condition number"}
    if (norm == "1"){Ylab = "condition number under 1-norm"}
    if (digitLoss | rlDist){par(mar = c(5,4,4,5)+.1)}
    plot(log(lambdas), type = "l", condNR, axes = FALSE, col = "blue4",
         xlab = "ln(penalty value)", ylab = Ylab, main = Main)
    axis(2, ylim = c(0,max(condNR)), col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    minor.tick(nx = 10, ny = 0, tick.ratio = .4)
    par(xpd = FALSE)
    if (digitLoss){
      dLoss <- floor(log10(condNR))
      par(new = TRUE)
      plot(log(lambdas), dLoss, axes = FALSE, type = "l", col = "green3",
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      axis(4, col = "black", lwd = 1)
      mtext("Loss in digits of accuracy", side = 4, line = 3)
      legend("top", col=c("blue4","green3"), lty = 1, legend =
               c("Condition number", "floor(log10(Condition number))"), cex = .8)
    }
    if (rlDist){
      RlDist <- 1/condNR
      par(new = TRUE)
      plot(log(lambdas), RlDist, axes = FALSE, type = "l", col = "green3",
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      axis(4, col = "black", lwd = 1)
      mtext("relative distance to singular matrix", side = 4, line = 3)
      legend("top", col=c("blue4","green3"), lty = 1, legend =
               c("Condition number", "Relative distance"), cex = .8)
    }
    if (vertical){
      if (missing(value)){
        stop("Need to specify input (value)")
      } else if (class(value) != "numeric"){
        stop("Input (value) is of wrong class")
      } else if (length(value) != 1){
        stop("Input (value) must be a scalar")
      } else if (value <= 0){
        stop("Input (value) must be positive")
      } else {
        abline(v = log(value), col = "red")
      }
    }

    # Possible output
    if (nOutput){
      return(list(lambdas = lambdas, conditionNumbers = condNR))
    }
  }
}



ridgeS <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Function that calculates Ridge estimators of a precision matrix
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a precision
  #             matrix
  # - ArchI   > Archetypal I ridge estimator of a precision matrix
  # - ArchII  > Archetypal II ridge estimator of a precision matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II est.
  # - When target is not the null-matrix it is expected to be p.d. for the
  #   vWP type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  # Dependencies
  # require("base")
  # require("expm")

  # Deprecation warning
  warning("This function is deprecated. Please use ridgeP instead.")

  if (!isSymmetric(S)) {
    stop("S should be a symmetric matrix")
  }
  else if (lambda <= 0) {
    stop("lambda should be positive")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else{
    # Calculate Ridge estimator
    # Alternative estimator
    if (type == "Alt"){
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (!all(target == 0) &
                 any(eigen(target, symmetric = TRUE,
                           only.values = TRUE)$values <= 0)){
        stop("When target is not a null-matrix it should be p.d. for this ",
             "type of ridge estimator")
      } else if (all(target == 0)){
        Spectral  <- eigen(S, symmetric = TRUE)
        Eigshrink <- .eigShrink(Spectral$values, lambda)
        P_Alt     <- solve(Spectral$vectors %*% diag(Eigshrink) %*%
                             t(Spectral$vectors))
        colnames(P_Alt) = rownames(P_Alt) <- colnames(S)
        return(P_Alt)
      } else if (all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
        varPhi    <- unique(diag(target))
        Spectral  <- eigen(S, symmetric = TRUE)
        Eigshrink <- .eigShrink(Spectral$values, lambda, const = varPhi)
        P_Alt     <- solve(Spectral$vectors %*% diag(Eigshrink) %*%
                             t(Spectral$vectors))
        colnames(P_Alt) = rownames(P_Alt) <- colnames(S)
        return(P_Alt)
      } else {
        E     <- (S - lambda * target)
        P_Alt <- solve(E/2 + sqrtm((E %*% E)/4 + lambda * diag(nrow(S))))
        return(P_Alt)
      }
    }

    # Archetypal I
    if (type == "ArchI"){
      if (lambda > 1){
        stop("lambda should be in (0,1] for this type of Ridge estimator")
      } else if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (any(eigen(target, symmetric = TRUE,
                           only.values = TRUE)$values <= 0)){
        stop("Target should always be p.d. for this type of ridge estimator")
      } else {
        P_ArchI <- solve((1-lambda) * S + lambda * solve(target))
        return(P_ArchI)
      }
    }

    # Archetypal II
    if (type == "ArchII"){
      P_ArchII <- solve(S + lambda * diag(nrow(S)))
      return(P_ArchII)
    }
  }
}



optPenalty.LOOCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt",
                             cor = FALSE, target = default.target(covML(Y)),
                             output = "light", graph = TRUE, verbose = TRUE) {
  ##############################################################################
  # - Function that selects the optimal penalty parameter by leave-one-out
  #   cross-validation
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax   > maximum value penalty parameter (dependent on 'type')
  # - step        > determines the coarseness in searching the grid
  #                 [lambdaMin, lambdaMax]
  # - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - cor         > logical indicating if evaluation of the LOOCV score should be
  #                 performed on the correlation matrix
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - output      > must be one of {"all", "light"}, default = "light"
  # - graph       > Optional argument for visualization optimal penalty
  #                 selection, default = TRUE
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("graphics")
  # require("sfsmisc")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (class(cor) != "logical"){
    stop("Input (cor) is of wrong class")
  }
  else if (!(output %in% c("all", "light"))){
    stop("Input (output) should be one of {'all', 'light'}")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else {
    # Set preliminaries
    LLs     <- numeric()
    lambdas <- lseq(lambdaMin, lambdaMax, length = step)

    # Calculate CV scores
    if (verbose) {
      cat("Calculating cross-validated negative log-likelihoods...\n")
    }
    for (k in 1:length(lambdas)){
      slh <- numeric()
      for (i in 1:nrow(Y)){
        S   <- covML(Y[-i,], cor = cor)
        slh <- c(slh, .LL(t(Y[i,,drop = F]) %*% Y[i,,drop = F],
                          ridgeP(S, lambdas[k], type = type, target = target)))
      }

      LLs <- c(LLs, mean(slh))
      if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
    }

    # Visualization
    optLambda <- min(lambdas[which(LLs == min(LLs))])
    if (graph){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
      plot(log(lambdas), type = "l", LLs, axes = FALSE,
           xlab = "ln(penalty value)", ylab = "LOOCV neg. log-likelihood",
           main = Main)
      axis(2, ylim = c(min(LLs),max(LLs)), col = "black", lwd = 1)
      axis(1, col = "black", lwd = 1)
      par(xpd = FALSE)
      abline(h = min(LLs), v = log(optLambda), col = "red")
      legend("topright",
             legend = c(paste("min. LOOCV neg. LL: ", round(min(LLs),3),sep=""),
                        paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
    }

    # Return
    S <- covML(Y, cor = cor)
    if (output == "all"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target),
                  lambdas = lambdas, LLs = LLs))
    }
    if (output == "light"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target)))
    }
  }
}



optPenalty.LOOCVauto <- function(Y, lambdaMin, lambdaMax,
                                 lambdaInit = (lambdaMin + lambdaMax)/2,
                                 cor = FALSE, target = default.target(covML(Y)),
                                 type = "Alt") {
  ##############################################################################
  # - Function that determines the optimal value of the penalty parameter by
  #   application of the Brent algorithm to the (leave-one-out) cross-validated
  #   log-likelihood
  # - Y          > (raw) Data matrix, variables in columns
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - cor        > logical indicating if evaluation of the LOOCV score should be
  #                performed on the correlation matrix
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  # input checks
  if (!is.matrix(Y)){
    stop("Input (Y) if of wrong class")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input matrix (Y) should not contain missings")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be strictly positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(lambdaInit) != "numeric"){
    stop("Input (lambdaInit) is of wrong class")
  }
  else if (length(lambdaInit) != 1){
    stop("Input (lambdaInit) must be a scalar")
  }
  else if (lambdaInit <= lambdaMin){
    stop("Input (lambdaInit) must be larger than input (lambdaMin)")
  }
  else if (lambdaMax <= lambdaInit){
    stop("Input (lambdaInit) must be smaller than input (lambdaMax)")
  }
  else if (class(cor) != "logical"){
    stop("Input (cor) is of wrong class")
  }
  else {
    # determine optimal value of ridge penalty parameter
    optLambda <- optim(lambdaInit, .cvl, method = "Brent", lower = lambdaMin,
                       upper = lambdaMax, Y = Y, cor = cor, target = target,
                       type = type)$par

    # Return
    return(list(optLambda = optLambda,
                optPrec = ridgeP(covML(Y, cor = cor), optLambda,
                                 type = type, target = target)))
  }
}
