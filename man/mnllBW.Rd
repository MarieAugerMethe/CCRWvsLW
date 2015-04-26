\name{mnllBW}
\alias{mnllBW}
\alias{nllBW}
\docType{package}
\title{Maximum likelihood estimation of the BW}
\description{Minimize the negative log-likelihood for the BW.}
\usage{
mnllBW(SL, TA, SLmin=min(SL))
nllBW(SL,TA,lambda,parF=list('SLmin'=min(SL)))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length, default min(SL)}
  \item{lambda}{one numeric value for lambda}
  \item{parF}{list for fix parameter values, see default}
}
\details{The mnllBW function uses an analytical solution to maximize the likelihood to estimate the parameters.

The nllBW function evaluates the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. arXiv:1406.4355 found at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmBW(500,0.1,1)
formPath <- movFormat(simPath)
mnllBW(formPath$SL, formPath$TA, formPath$SLmin)
}
