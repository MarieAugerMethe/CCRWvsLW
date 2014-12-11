\name{simmBW}
\alias{simmBW}
\title{Simulate a Brownian Walk (BW)}
\description{Simulate a Brownian Walk (BW). This model represent the movement of an animal which has no directional persistence.}
\usage{
simmBW(n, l, a)
}
\arguments{
  \item{n}{one integer value for the sample size. Note that this sample size represent the number of final step lengths and turning angles wanted. The ltraj object returned will be longer because a minimum of 3 locations are required to calculate a relative turning angle}
  \item{l}{one numeric and positive value representing the lambda value for the exponential distribution of the step lengths}
  \item{a}{one numeric and positive value representing the minimum step length value}
}
\details{Maybe talk about Null model}
\value{SAY WHAT IT RETUNS (ltraj object!)}
\references{MY PAPER}