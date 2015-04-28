\name{emHMM}
\alias{emHMM}
\alias{EMHMM}
\title{EM-algorithm to fit a hidden Markov model representing the CCRW
}
\description{
emHMM finds the maximum likelihood estimate for the parameters of the CCRW by fitting the hidden Markov model through an Expectation Maximization (EM) algorithm.
}
\usage{
emHMM(SL, TA, missL, SLmin, lambda, gamm, delta=c(0.5,0.5), 
kapp, notMisLoc, maxiter=10000, tol=1e-5)
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{missL}{integer vector containing the number of time steps between two steps. If no missing location it will be 1.}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{lambda}{numeric vector of length 2 containing the starting value for the lambdas of the two behaviors}
  \item{gamm}{2x2 matrix containing the starting value for the transition probability matrix}
  \item{delta}{numeric vector value for the probability of starting in each of the two behaviors, default value c(0.5,0.5), which means that you have an equal chance of starting in each behavior}
  \item{kapp}{one numeric value representing the starting value for the kappa of the von Mises distribution describing the extensive search behavior}
  \item{notMisLoc}{integer vector containing the index of the locations that are not missing}
  \item{maxiter}{one integer value representing the maximum number of iterations the EM algorithm will go through. Default = 10000.}
  \item{tol}{double: value that indicates the maximum allowed difference between the parameters.}
  }
\details{Will return the parameter estimates and the minimum negative log likelihood.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}

For more information on the EM-algorithm please refer to Zucchini W. and I.L. MacDonald (2009) Hidden Markov Models for Time Series: An Introduction Using R. Chapman and Hall/CRC}
\examples{
simPath <- simmCCRW(500,0.9,0.9,0.1,0.01,5,1)
formPath <- movFormat(simPath)
emHMM(formPath$SL,formPath$TA, formPath$missL, formPath$SLmin, 
      lambda=c(0.1,0.1), gamm=matrix(c(0.8,0.2,0.2,0.8),nrow=2),
      kapp=10,notMisLoc=formPath$notMisLoc)
}
