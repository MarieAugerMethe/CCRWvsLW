\name{mnllTCRW}
\alias{mnllTCRW}
\docType{package}
\title{
Muck.
}
\description{
Muck. 
}
\usage{
mnllTCRW(SL, TA_C, TA_N, SLmin, SLmax,
  lambda=exp(nlm(nllTBW, log(1/mean(SL-SLmin)),SL=SL,TA_N=TA_N,
                 SLmin=SLmin,SLmax=SLmax)$estimate),
	kapp=mle.vonmises(TA_C,mu=circular(0))$kappa)
}
\arguments{
  \item{SL}{
    double: muck
  }
  \item{TA_C}{
    double: muck
  }
  \item{TA_N}{
    double: muck
  }
  \item{SLmin}{
    double: muck
  }
  \item{SLmax}{
    double: muck
  }
  \item{lambda}{
    double: muck
  }
  \item{kapp}{
    double: muck
  }
}
