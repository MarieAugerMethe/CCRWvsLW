\name{mnllCCRW}
\alias{mnllCCRW}
\title{Maximum likelihood estimation of the CCRW}
\description{Minimise the negative log-likelihood for the CCRW though the EM-algorithm.
}
\usage{
mnllCCRW(SL, TA, TA_C, missL, notMisLoc, SLmin, tol=5e-5)
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{TA_C}{circular object containing the turning angles}
  \item{missL}{integer vector containing the number of time step between two steps. If no missing location it will be 1.}
  \item{notMisLoc}{integer vector containing the index of the locations that are not missing}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{tol}{double: value that indicates the maximum allowed difference between the parameters, default 5e-5}
}
\details{This function use the function emHMM to minimize the negative log-likelihood. Note that this rquire setting starting paramater values. To limit the chances to that the estimates is the results of a local peak in the likelihood, the function mnllCCRW applies the EM-algorithm function emHMM to set of different starting values combinations and choose the overall lowest negative log-likelihood estimate. To explore further starting parameter values use directly emHMM function.}
\seealso{\code{\link{emHMM}}}
