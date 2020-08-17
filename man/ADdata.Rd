\name{ADdata}
\alias{ADdata}
\alias{ADmetabolites}
\alias{sampleInfo}
\alias{variableInfo}
\docType{data}
\title{
R-objects related to metabolomics data on patients with Alzheimer's Disease
}
\description{
\code{ADdata} contains 3 objects related to metabolomics data on patients with Alzheimer's Disease (AD).

\code{ADmetabolites} is a \code{matrix} containing metabolic expressions of 230 metabolites (rows) on 127 samples (columns).

\code{sampleInfo} is a \code{data.frame} containing information on the samples.
Information pertains to diagnosis: AD class 1 or AD class 2.

\code{variableInfo} is a \code{data.frame} containing information on the metabolic features.
Information pertains to compound families: Amines, organic acids, lipids, and oxidative stress compounds.
}
\usage{data(ADdata)}
\details{
See description.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\source{
de Leeuw, F., Peeters, C.F.W., Kester, M.I., Harms, A.C., Struys, E., Hankemeijer, T., van Vlijmen, H.W.T., van Duijn, C.M., Scheltens, P., Demirkan, A., van de Wiel, M.A., van der Flier, W.M., and Teunissen, C.E. (2017). Blood-based metabolic signatures in Alzheimer's Disease. Alzheimer's & Dementia: Diagnosis, Assessment & Disease Monitoring, 8: 196-207.
}
\examples{
data(ADdata)

## Look at sample information
sampleInfo

## Look at feature information
variableInfo
}
\keyword{datasets}
