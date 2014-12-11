\name{simmCRW}
\alias{simmCRW}
\title{Simulate a Correlated Random Walk (CRW)}
\description{Simulate a Correlated Random Walk (CRW). This model represent the movement of an animal which has directional persistence.}
\usage{
simmCRW(n, l, k, a)
}
\arguments{
  \item{n}{one integer value for the sample size. Note that this sample size represent the number of final step lengths and turning angles wanted. The ltraj object returned will be longer because a minimum of 3 locations are required to calculate a relative turning angle}
  \item{l}{one numeric and positive value representing the lambda value for the exponential distribution of the step lengths}
  \item{k}{one numeric and positive value representing the kappa value of von Mises distribution of the turning angles}
  \item{a}{one numeric and positive value representing the minimum step length value}
}
\details{Maybe talk about Null model}
\value{SAY WHAT IT RETUNS (ltraj object!)}
\references{MY PAPER}