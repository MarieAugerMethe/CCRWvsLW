\name{movFormat}
\alias{movFormat}
\docType{package}
\title{Transform a one path ltraj object into the measures used to do the movement analysis.
}
\description{
This function takes a ltraj object and calculates the step length and turning angle as well as the missing locations.
}
\usage{
movFormat(movltraj, TAc = 0)
}
\arguments{
  \item{movltraj}{ltraj object}
  \item{TAc}{Threshold angle used in the local turn method}
}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
