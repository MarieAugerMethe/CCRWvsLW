\name{movLikelihoods}
\alias{movLikelihoods}
\docType{package}
\title{Fit seven random walk models to a movement path}
\description{Minimize the negative log-likelihood for the BW.}
\usage{
movLikelihoods(movltraj, graph = TRUE, PRdetails = FALSE, TAc = 0, conts=TRUE, dn=FALSE, ww=FALSE, hs=FALSE)
}
\arguments{
  \item{movltraj}{one ltraj object containing the movement path}
  \item{graph}{one logical value that states whether the models should be plotted}
  \item{PRdetails}{one logical value that states whether the details of the pseudoresiduals should be returned}
  \item{TAc}{one numeric value for the threshold angle used in the local turn method}
  \item{conts}{logical value stating whether to contrain mu between 1 and 3 (i.e. the values relevant for the Levy walk searching strategy)}
  \item{dn}{logical value stating whether to direct numerical maximization of the log likelihood is used rather than the EM-algorithm for the CCRW. see details}
  \item{ww}{logical value stating whether to investigate the CCRW with weibull distributions for SL and wrapped cauchy for TA.}
  \item{hs}{logical value stating whether to investigate the CCRW that is a hidden semi Markov model (HSMM) with weibull distributions for SL and wrapped cauchy for TA.}
}
\details{movLikelihoods fits a CCRW, TLW, LW, BW, CRW, TBW, TCRW to a movement path.

Note by default the EM-algorithm is used for the CCRW. Because it is coded in C++ it is much faster than the numerical maximization. If dn=TRUE, the numerical maximization is used, which has the advantage to have one less parameter to estimate: dI.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015-In Press) Differentiating the Levy walk from a composite correlated random walk. Methods in Ecology and Evolution. Preprint available at \url{http://arxiv.org/abs/1406.4355}
}
\examples{
simPath <- simmCCRW(500,0.9,0.9,0.1,0.01,5,1)
res <- movLikelihoods(simPath)
res$mleMov
AICcRes <- unlist(res$mle)[grep("AICc", names(unlist(res$mle)))]
AICcRes-min(AICcRes)
res$CI
res$pseudoRes
}
