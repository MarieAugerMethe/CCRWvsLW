\name{mnllTBW}
\alias{mnllTBW}
\alias{nllTBW}
\docType{package}
\title{Maximum likelihood estimation of the TBW}
\description{Minimize the negative log-likelihood for the TBW.}
\usage{
mnllTBW(SL, TA, SLmin=min(SL), SLmax=max(SL), lambda_BW=(1/mean(SL-SLmin)))
nllTBW(SL,TA,x,parF=list('SLmin'=SLmin,'SLmax'=SLmax))
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{SLmin}{one numeric value representing the minimum step length, default min(SL)}
  \item{SLmax}{one numeric value representing the maximum step length, default max(SL)}
  \item{lambda_BW}{one numeric value for the starting value of lambda, the default is the analytical solution BW}
  \item{x}{one numeric value for log lambda}
  \item{parF}{list for fix parameter values, see default}
}
\details{The mnllTBW function uses the value the analytical solution for the BW as the starting value of the numerical minimser to get lambda for the TBW. It returns the parameter values, the minimum negative log likelihood, AIC and AICc values. 

The nllTBW function evaluates the negative log likelihood value.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmBW(500,0.1,1)
formPath <- movFormat(simPath)
mnllTBW(formPath$SL, formPath$TA)
}