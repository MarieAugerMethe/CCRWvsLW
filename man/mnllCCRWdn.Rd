\name{mnllCCRWdn}
\alias{mnllCCRWdn}
\title{Numerical maximum likelihood estimation of the CCRW}
\description{Minimize the negative log-likelihood for the CCRW though numerical maximization.
}
\usage{
mnllCCRWdn(SL, TA, TA_C, missL, SLmin)
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{TA_C}{circular object containing the turning angles}
  \item{missL}{integer vector containing the number of time step between two steps. If no missing location it will be 1.}
  \item{SLmin}{one numeric value representing the minimum step length}
}
\details{This function uses the optim to numerically minimize the negative log-likelihood. Note that optim requires setting starting parameter values. To decrease the chances of having results associated with a local peak in the likelihood, the function mnllCCRW applies use set of different starting values combinations and choose the overall lowest negative log-likelihood estimate. To explore the effect of starting values, use directly the nllCCRWdn function.}
\seealso{\code{\link{mnllCCRW}}}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCCRW(500,0.9,0.9,0.1,0.01,5,1)
formPath <- movFormat(simPath)
mnllCCRWdn(formPath$SL, formPath$TA, formPath$TA_C, 
      formPath$missL, formPath$SLmin)
}
