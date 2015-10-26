\name{mnllTLW}
\alias{mnllTLW}
\alias{nllTLW}
\docType{package}
\title{Maximum likelihood estimation of the TLW}
\description{Minimize the negative log-likelihood for the TLW.}
\usage{
mnllTLW(SL, TA, SLmin=min(SL), SLmax=max(SL), conts=TRUE)
nllTLW(SL,TA,mu,parF=list(SLmin=min(SL),SLmax=max(SL)))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length, default is the min(SL)}
  \item{SLmax}{one numeric value representing the maximum step length, default is the max(SL)}
  \item{mu}{one numeric value for parameter}
  \item{parF}{list for fix parameter values, see default}
  \item{conts}{logical value stating whether to contrain mu between 1 and 3 (i.e. the values relevant for the Levy walk searching strategy)}
}
\details{The mnllTLW function numerically minimeze the negative log likelihood to estimate the parameter mu. Note that the parameter value is constrained between 1 and 3.

The nllTLW function evaluates the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmTLW(500,2,1,1000)
formPath <- movFormat(simPath)
mnllTLW(formPath$SL, formPath$TA)
}