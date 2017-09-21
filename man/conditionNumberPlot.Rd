\name{conditionNumberPlot}
\alias{conditionNumberPlot}
\title{
Visualize the spectral condition number against the regularization parameter
}
\description{
This function is now deprecated. Please use \code{CNplot} instead.
}
\usage{
conditionNumberPlot(S, lambdaMin, lambdaMax, step, type = "Alt",
                    target = default.target(S), norm = "2",
                    digitLoss = FALSE, rlDist = FALSE,
                    vertical = FALSE, value, main = TRUE,
                    nOutput = FALSE, verbose = TRUE)
}
\arguments{
  \item{S}{
Sample covariance \code{matrix}.
}
  \item{lambdaMin}{
A \code{numeric} giving the minimum value for the penalty parameter.
}
  \item{lambdaMax}{
A \code{numeric} giving the maximum value for the penalty parameter.
}
  \item{step}{
An \code{integer} determining the number of steps in moving through the grid [\code{lambdaMin}, \code{lambdaMax}].
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{norm}{
A \code{character} indicating the norm under which the condition number is to be calculated/estimated. Must be one of: "1", "2".
}
  \item{digitLoss}{
A \code{logical} indicating if the approximate loss in digits of accuracy should also be visualized in the output graph.
}
  \item{rlDist}{
A \code{logical} indicating if the relative distance to the set of singular matrices should also be visualized in the output graph.
}
  \item{vertical}{
A \code{logical} indicating if output graph should come with a vertical line at a pre-specified value for the penalty parameter.
}
  \item{value}{
A \code{numeric} indicating a pre-specified value for the penalty parameter.
}
  \item{main}{
A \code{logical} indicating if output graph should contain type of estimator as main title.
}
  \item{nOutput}{
A \code{logical} indicating if numeric output should be returned.
}
  \item{verbose}{
A \code{logical} indicating if information on progress should be printed on screen.
}
}
\details{
See \code{CNplot}.
}
\value{
The function returns a graph. If \code{nOutput = TRUE} the function also returns an object of class \code{list}:
\item{lambdas}{A \code{numeric} vector representing all values of the penalty parameter for which the condition number was calculated.}
\item{conditionNumbers}{A \code{numeric} vector containing the condition number for each value of the penalty parameter given in \code{lambdas}.}
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\seealso{
\code{\link{CNplot}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Assess spectral condition number across grid of penalty parameter
conditionNumberPlot(Cx, lambdaMin = .0001, lambdaMax = 50, step = 1000)
}
