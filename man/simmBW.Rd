\name{simmBW}
\alias{simmBW}
\title{Simulate a Brownian Walk (BW)}
\description{Simulate a Brownian Walk (BW). This model represent the movement of an animal which has no directional persistence.}
\usage{
simmBW(n, l, a)
}
\arguments{
  \item{n}{one integer value for the sample size. Note that this sample size represent the number of final step lengths and turning angles wanted (when TAc = 0). The ltraj object returned will be longer because a minimum of 3 locations are required to calculate a relative turning angle}
  \item{l}{one numeric and positive value representing the lambda value for the exponential distribution of the step lengths}
  \item{a}{one numeric and positive value representing the minimum step length value}
}
\details{Simulates a BW and return a ltraj object}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
