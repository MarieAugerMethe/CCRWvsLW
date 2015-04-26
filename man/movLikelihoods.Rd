\name{movLikelihoods}
\alias{movLikelihoods}
\docType{package}
\title{Fit seven random walk models to a movement path}
\description{Minimize the negative log-likelihood for the BW.}
\usage{
movLikelihoods(movltraj, graph = TRUE, PRdetails = FALSE, TAc = 0)
}
\arguments{
  \item{movltraj}{one ltraj object containing the movement path}
  \item{graph}{one logical value that states whether the models should be plotted}
  \item{PRdetails}{one logical value that states whether the details of the pseudoresiduals should be returned}
  \item{TAc}{one numeric value for the threshold angle used in the local turn method}
}
\details{movLikelihoods fits a CCRW, TLW, LW, BW, CRW, TBW, TCRW to a movement path.}
\references{
Please refer to Auger-Methe, M., A.E. Derocher, M.J. Plank, E.A. Codling, M.A. Lewis (2015) Differentiating the Levy walk from a composite correlated random walk. arXiv:1406.4355 found at \url{http://arxiv.org/abs/1406.4355}
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
