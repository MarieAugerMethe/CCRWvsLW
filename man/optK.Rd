\name{optK}
\alias{optK}
\docType{package}
\title{
Minimizing numerically the MLE to estimate kappa.
}
\description{
Even with the CDLL the kappa value cannot be minimized analytically. This function minimize the part that needs to be minimized numerically. 
}
\usage{
optK(t2)	
}
\arguments{
  \item{t2}{double: part of the derivative, see details}
}
\details{The get maximum likelihood estimate for kappa using the CDLL, we need to solve the following equation: 0 = I_1(kappa)/I_0(kappa) - sum (E(u_j(t))*cos(x_t))/sum(E(u_j(t))), where I_1 and I_0 are the modified Bessel functions. This equation is hard to solve analytically and optK uses the function optimize to numerically minimize: abs((besselI(k, 1)/besselI(k,0)) - t2), where t2 is: sum (E(u_j(t))*cos(x_t))/sum(E(u_j(t))). This is only done for the extensive behavior (j=2). The kappa of the intensive behavior (j=1) is fixed to 0 (uniform distribution).
}