\name{simmTLW}
\alias{simmTLW}
\title{Simulate a truncated Levy Walk (TLW)}
\description{Simulate a truncated Levy Walk (TLW). This model represent the movement of the Levy search strategy.}
\usage{
simmTLW(n, mu, a, b)
}
\arguments{
  \item{n}{one integer value for the sample size. Note that this sample size represent the number of final step lengths and turning angles wanted (when TAc = 0). The ltraj object returned will be longer because a minimum of 3 locations are required to calculate a relative turning angle}
  \item{mu}{one numeric (1<mu<=3) representing the mu value for the Pareto distribution of the step lengths}
  \item{a}{one numeric and positive value representing the minimum step length value}
  \item{b}{one numeric and positive value representing the maximum step length value, i.e., the truncation point}
}
\details{Maybe talk about Null model}
\value{SAY WHAT IT RETUNS (ltraj object!)}
\references{MY PAPER}
\seealso{\code{\link{simmLW}}}