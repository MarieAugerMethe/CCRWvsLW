\name{mnllCCRW}
\alias{mnllCCRW}
\title{Maximum likelihood estimation of the CCRW}
\description{Minimize the negative log-likelihood for the CCRW though the EM-algorithm.
}
\usage{
mnllCCRW(SL, TA, TA_C, missL, notMisLoc, SLmin, tol=5e-5, conts=TRUE)
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{TA_C}{circular object containing the turning angles}
  \item{missL}{integer vector containing the number of time step between two steps. If no missing location it will be 1.}
  \item{notMisLoc}{integer vector containing the index of the locations that are not missing}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{tol}{double: value that indicates the maximum allowed difference between the parameters, default 5e-5}
  \item{conts}{logical value stating whether to contrain mu between 1 and 3 (i.e. the values relevant for the Levy walk searching strategy)}
}
\details{This function uses the function emHMM to minimize the negative log-likelihood. Note that emHMM requires setting starting parameter values. To decrease the chances of having results associated with a local peak in the likelihood, the function mnllCCRW applies the EM-algorithm function emHMM to set of different starting values combinations and choose the overall lowest negative log-likelihood estimate. To explore the effect of starting values, use directly emHMM function.}
\seealso{\code{\link{emHMM}}}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCCRW(500,0.9,0.9,0.1,0.01,5,1)
formPath <- movFormat(simPath)
mnllCCRW(formPath$SL, formPath$TA, formPath$TA_C, 
      formPath$missL, formPath$notMisLoc, formPath$SLmin)
}
