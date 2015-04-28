\name{mnllCRW}
\alias{mnllCRW}
\alias{nllCRW}
\docType{package}
\title{Maximum likelihood estimation of the CRW}
\description{Minimize the negative log-likelihood for the CRW.}
\usage{
mnllCRW(SL, TA_C, TA, SLmin, lambda=(1/mean(SL-SLmin)))
nllCRW(SL,TA,kapp,lambda=parF$lambda,parF=list('SLmin'=min(SL)))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{TA_C}{circular object containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{lambda}{one numeric value for lambda, the default is the analytical solution}
  \item{kapp}{one numerci value for parameter}
  \item{parF}{list for fix parameter values, see default}
}
\details{The mnllCRW function minimizes the negative log likelihood for the CRW. It uses the value the analytical solution that is the input for lambda and gets the MLE for kappa using the mle.vonmises function of the package circular.

The nllCRW function evaluates the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCRW(500,0.1,5,1)
formPath <- movFormat(simPath)
mnllCRW(formPath$SL, formPath$TA_C, formPath$TA, formPath$SLmin)
}