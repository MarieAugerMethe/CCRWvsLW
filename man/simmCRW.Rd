\name{simmCRW}
\alias{simmCRW}
\title{Simulate a Correlated Random Walk (CRW)}
\description{Simulate a Correlated Random Walk (CRW). This model represent the movement of an animal which has directional persistence.}
\usage{
simmCRW(n, l, k, a)
}
\arguments{
  \item{n}{one integer value for the sample size. Note that this sample size represent the number of final step lengths and turning angles wanted (when TAc = 0). The ltraj object returned will be longer because a minimum of 3 locations are required to calculate a relative turning angle}
  \item{l}{one numeric and positive value representing the lambda value for the exponential distribution of the step lengths}
  \item{k}{one numeric and positive value representing the kappa value of von Mises distribution of the turning angles}
  \item{a}{one numeric and positive value representing the minimum step length value}
}
\details{Simulates a CRW and return a ltraj object}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. arXiv:1406.4355 found at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCRW(500,0.1,15,1)
simPath
plot(simPath)
}