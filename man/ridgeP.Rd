\name{ridgeP}
\alias{ridgeP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Same as ridgeS
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ridgeP(S, lambda, type = "Alt", target = default.target(S))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (S, lambda, type = "Alt", target = default.target(S))
{
    if (!isSymmetric(S)) {
        stop("S should be a symmetric matrix")
    }
    else if (lambda <= 0) {
        stop("lambda should be positive")
    }
    else if (!(type \%in\% c("Alt", "ArchI", "ArchII"))) {
        stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
    }
    else {
        if (type == "Alt") {
            if (!isSymmetric(target)) {
                stop("Shrinkage target should be symmetric")
            }
            else if (dim(target)[1] != dim(S)[1]) {
                stop("S and target should be of the same dimension")
            }
            else {
                P_Alt <- .armaRidgeP(S, target, lambda)
            }
            dimnames(P_Alt) <- dimnames(S)
            return(P_Alt)
        }
        if (type == "ArchI") {
            if (lambda > 1) {
                stop("lambda should be in (0,1] for this type of Ridge estimator")
            }
            else if (!isSymmetric(target)) {
                stop("Shrinkage target should be symmetric")
            }
            else if (dim(target)[1] != dim(S)[1]) {
                stop("S and target should be of the same dimension")
            }
            else if (any(eigen(target, symmetric = TRUE, only.values = TRUE)$values <=
                0)) {
                stop("Target should always be p.d. for this type of ridge estimator")
            }
            else {
                P_ArchI <- solve((1 - lambda) * S + lambda *
                  solve(target))
                return(P_ArchI)
            }
        }
        if (type == "ArchII") {
            P_ArchII <- solve(S + lambda * diag(nrow(S)))
            return(P_ArchII)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
