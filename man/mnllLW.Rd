\name{mnllLW}
\alias{mnllLW}
\alias{nllLW}
\docType{package}
\title{Maximum likelihood estimation of the LW}
\description{Minimize the negative log-likelihood for the LW.}
\usage{
mnllLW(SL, TA, SLmin = min(SL))
nllLW(SL,TA,mu,parF=list('SLmin'=min(SL)))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths, default min(SL)}
  \item{TA}{numeric vector containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{mu}{one numeric value for the parameter}
  \item{parF}{list for fix parameter values, see default}
}
\details{The mnllLW function uses an analytical solution to maximize the likelihood to estimate the parameter mu. Note that the parameter value is restrained between 1 and 3, if the value returned by the analytical formula is greater outside these bounds a constrained optimizer is used.

The nllLW function evaluate the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. arXiv:1406.4355 found at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmLW(500,2,1)
formPath <- movFormat(simPath)
mnllLW(formPath$SL, formPath$TA)
}
