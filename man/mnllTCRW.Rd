\name{mnllTCRW}
\alias{mnllTCRW}
\alias{nllTCRW}
\docType{package}
\title{Maximum likelihood estimation of the TCRW}
\description{Minimize the negative log-likelihood for the TCRW.}
\usage{
mnllTCRW(SL, TA_C, TA, SLmin=min(SL), SLmax=max(SL),
    lambda = exp(nlm(nllTBW, log(1/mean(SL - SLmin)), 
    SL = SL, TA = TA, SLmin = SLmin, SLmax = SLmax)$estimate), 
    kapp = mle.vonmises(TA_C, mu = circular(0))$kappa)
nllTCRW(SL, TA, lambda, kapp, SLmin=min(SL), SLmax=max(SL))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA_C}{circular object containing the turning angles}
  \item{TA}{numeric vector containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{SLmax}{one numeric value representing the maximum step length}
  \item{lambda}{one numeric value for lambda value, the default MLE for the TBW}
  \item{kapp}{one numeric value for kappa value, the default MLE for the CRW}
}
\details{The mnllTCRW function minimize the negative log likelihood using the default formula in the input to get the estimates. It returns parameter estimates, the negative log likelihood, AIC, and AICc.

The nllTCRW function evaluate the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. arXiv:1406.4355 found at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCRW(500,0.1,5,1)
formPath <- movFormat(simPath)
mnllCRW(formPath$SL, formPath$TA_C, formPath$TA, formPath$SLmin)
}
