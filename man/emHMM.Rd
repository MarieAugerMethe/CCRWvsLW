\name{emHMM}
\alias{emHMM}
\alias{EMHMM}
\alias{toMinK}
\alias{optK}
\title{EM-algorithm to fit the hiden Markov model for the CCRW
}
\description{
emHMM finds the maximum likelihood estimate for the parameters of the CCRW by fitting the hiden Markov model through an Expectation Maximazation (EM) algorithm.
}
\usage{
emHMM(SL, TA, missL, SLmin, lambda, gamm, delta=c(0.5,0.5), kapp, notMisLoc, maxiter=10000, tol=1e-5)
}
\arguments{
  \item{SL}{numeric vector containing the step lengths}
  \item{TA}{numeric vector containing the turning angles}
  \item{missL}{integer vector containing the number of time step between two steps. If no missing location it will be 1.}
  \item{SLmin}{one numeric value representing the minimum step length}
  \item{lambda}{numeric vector of length 2 containging the starting value for the lambdas of the two behaviours}
  \item{gamm}{2x2 matrix containing the starting value for the transition probability matrix}
  \item{delta}{numeric vector value for the probability of starting in each of the two behaviours, default value c(0.5,0.5), which means that you have an equal chance of starting in each behavior}
  \item{kapp}{one numeric value representing the starting value for the kappa of the von Mises distribution describing the extensive search behavior}
  \item{notMisLoc}{integer vector containing the index of the locations that are not missing}
  \item{maxiter}{one integer value representing the maximum number of iterations the EM algorithm will go through. Default = 10000.}
  \item{tol}{double: value that indicates the maximum allowed difference between the parameters.}
  }
\details{MAYBE TALK ABOUT EMHMM the cpp function, optK, and toMinK}
\references{ADD CITATION TO ZUC AND MAC AND MY PAPER}