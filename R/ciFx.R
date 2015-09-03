#######################################
# Created By: Marie Auger-Methe
# Date created: November 21, 2011
# Updated: April 3, 2013

# This script contains the functions for computing the quadratic approximation CI (using the Hessian).
# The fx for the quad. app. CI CI.Hessian() can be used directly
# to calculate the CI for all parameters estimated through mle of a model.
# However, when the model has only one parameter
# or independet parameter one can use directly the second derivative
# instead of the Hessian.
# In cases where there is an analytical solution for the CI
# we use the analitycal solution.

#######################################
CI.Hessian <- function(SL,TA, parMLE, trans.par, M, parF, NLL){

	# Create matrix for Hessian C.I. results
	parCI <- matrix(NA,nrow=length(parMLE),ncol=2)
	colnames(parCI) <- c("L95CI", "U95CI")
	rownames(parCI) <- names(parMLE)

	# Get a numerial estimate of the Hessian from the optimisation routine
	# but this is for the working parameter
	# the nlm often crashes for CCRW so can't use it for hessian
	H <- hessian(NLL, trans.par(parMLE), SL=SL,TA=TA,parF=parF)

	# Transform for real parameter
	# From Zucch & MacD 2009
	# G^{-1} = M' %*% H^{-1} %*% M
	# Where H is the hessian of the working parameter,
	# G is the hessian of the true parameters,
	# M is define by m_ij = d theta_j / d phi_i
	# (where theta is true parameter and phi are working parameter)

	# although I feel like there is an easier solution
	# For X: d plogis(phi) /d phi
	# exp(phi)/(1+exp(phi))^2
	# For lambda_F, lambda_T, kappa_T:
	# d exp(phi)/d phi
	# exp(phi): i.e the true parameter
	
	# I'm not sure why I get NaN for det(H)
	# In happens in the case of kappa = 0
	if(abs(det(H)) < 1e-300 | is.nan(det(H))){
		warning('The determinant of the hessian is 0
			or computationally abs(det(H))<1e-300.
			Or the determinant is not a number.
			Maybe link to unbouded likelihood problems(?)')
	}else{
		# G_inv is the Variance-Covariance matrix
		# See Bolker p. 263 & Zucc. & MacD. 2009 Ch3
		G_inv <- t(M) %*% solve(H, tol=1e-300) %*% M

		int <- qnorm(0.975)*sqrt(diag(G_inv))
	
		parCI[,1] <- parMLE - int
		parCI[,2] <- parMLE + int
	}
	return(parCI)	

}


#######################################
# CCRW - for EM-algorithm

ciCCRW <- function(SL,TA,SLmin,missL,notMisLoc,mleM){
  # Table for the CI
  CI <- matrix(NA, nrow=5, ncol=3)
  rownames(CI) <- names(mleM[1:5])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimates
  CI[,1] <- mleM[1:5]
  
	parF <- list('SLmin'=SLmin)
	trans.par <- function(x){
		x[1] <- qlogis(x[1])
		x[2] <- qlogis(x[2])
		x[3] <- log(x[3] - .Machine$double.xmin)
		x[4] <- log(x[4] - .Machine$double.xmin)
		x[5] <- log(x[5])
		return(x)
	}

	# According to Zucchini and MacDonald (2009)
	M <- diag(c(
	  CI[1,1]*(1-CI[1,1]),
	  CI[2,1]*(1-CI[2,1]),
	  CI[3,1], CI[4,1], CI[5,1]))

	# Note that although I put places for the deltas
	# I won't get any CIs




  if(any(is.na(mleM[1:7]))){
	  warning("The EM-algorithm gives NA values for the parameters")
	}else{
	  CI[,2:3] <- CI.Hessian(SL,TA, CI[,1], trans.par, M,
	                         parF=c(parF,mleM[6],mleM[7],list("missL"=missL)), nllCCRW)
	}
  
	return(CI)
}

#######################################
# CCRW - direct numerical minimization

ciCCRWdn <- function(SL,TA,SLmin,missL,mleM){
  # Table for the CI
  CI <- matrix(NA, nrow=5, ncol=3)
  rownames(CI) <- names(mleM[1:5])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimates
  CI[,1] <- mleM[1:5]
  
  parF <- list('SLmin'=SLmin,"missL"=missL)
  trans.par <- function(x){
    x[1] <- qlogis(x[1])
    x[2] <- qlogis(x[2])
    x[3] <- log(x[3] - .Machine$double.xmin)
    x[4] <- log(x[4] - .Machine$double.xmin)
    x[5] <- log(x[5])
    return(x)
  }
  
  # According to Zucchini and MacDonald (2009)
  M <- diag(c(
    CI[1,1]*(1-CI[1,1]),
    CI[2,1]*(1-CI[2,1]),
    CI[3,1], CI[4,1], CI[5,1]))
  
  if(any(is.na(mleM[1:7]))){
    warning("The optim return NA values for the parameters, so no CI calculated")
  }else{
    CI[,2:3] <- CI.Hessian(SL,TA, CI[,1], trans.par, M,
                           parF=parF, nllCCRWdn)
  }
  
  return(CI)
}

#######################################
# CCRW - direct numerical minimization weibull and wrapped cauchy

ciCCRWww <- function(SL,TA,missL,mleM){
  # Table for the CI
  CI <- matrix(NA, nrow=7, ncol=3)
  rownames(CI) <- names(mleM[1:7])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimates
  CI[,1] <- mleM[1:7]
  
  parF <- list("missL"=missL)
  trans.par <- function(x){
    x[1] <- qlogis(x[1])
    x[2] <- qlogis(x[2])
    x[3] <- log(x[3] - .Machine$double.xmin)
    x[4] <- log(x[4] - .Machine$double.xmin)
    x[5] <- log(x[5] - .Machine$double.xmin)
    x[6] <- log(x[6] - .Machine$double.xmin)
    x[7] <- qlogis(x[7])
    return(x)
  }
  
  M <- diag(c(
    CI[1,1]*(1-CI[1,1]),
    CI[2,1]*(1-CI[2,1]),
    CI[3,1], CI[4,1], CI[5,1],CI[6,1],CI[7,1]))
  
  if(any(is.na(mleM[1:7]))){
    warning("The optim return NA values for the parameters, so no CI calculated")
  }else{
    CI[,2:3] <- CI.Hessian(SL,TA, CI[,1], trans.par, M,
                           parF=parF, nllCCRWww)
  }
  
  return(CI)
}


#######################################
# CCRW - direct numerical minimization weibull and wrapped cauchy

ciHSMM <- function(SL,TA,missL,notMisLoc,mleM){
  # Table for the CI
  CI <- matrix(NA, nrow=9, ncol=3)
  rownames(CI) <- names(mleM[1:9])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimates
  CI[,1] <- mleM[1:9]
  
  parF <- list("missL"=missL, "notMisLoc"=notMisLoc, m=c(20,20))
  trans.par <- function(x){
    x[1:8] <-  log(x[1:8] - .Machine$double.xmin)
    x[9] <- qlogis(x[9])
    return(x)
  }
  
  M <- diag(CI[1:9,1])
  
  if(any(is.na(mleM[1:7]))){
    warning("The optim return NA values for the parameters, so no CI calculated")
  }else{
    CI[,2:3] <- CI.Hessian(SL,TA, CI[,1], trans.par, M,
                           parF=parF, nllHSMM)
  }
  
  return(CI)
}


#######################################
# LW

ciLW <- function(SL,TA,SLmin,mleM){
  # Table for the CI
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimate
  CI[,1] <- mleM[1]
  
	parF <- list('SLmin'=SLmin)
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# In the case of the LW
	# there is an analytical solution for the second derivative
	# D2 = -n/(mu-1)^2
	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	# 1/sqrt(-n/(mu-1)^2) = (mu-1)/sqrt(-n)
	int <- qnorm(0.975)*(CI[,1]-1)/sqrt(length(SL))
	CI[,2] <- CI[,1] - int
	CI[,3] <- CI[,1] + int

	return(CI)
}

#######################################
# TLW

ciTLW <- function(SL,TA,SLmin,SLmax,mleM){
  # Table for CI
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimate
  CI[,1] <- mleM[1]
  
	parF <- list('SLmin'=SLmin,'SLmax'=SLmax)
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	

	# In the case of the TLW
	# there is an analytical solution for the second derivative
	# But it is complcated.
	# So I'm using R to get the second derivative
	# I'm first deriving the negative log pdf in two parts
	# Part 1 doesn't include the SL 
	npdf_1 <- expression(-(log(mu-1) - log(a^(1-mu) - b^(1-mu))))
	D2_1 <- D(D(npdf_1,'mu'),'mu')
	# Part 2 includes the SL
	# but npdf_2 = 0 so I am not running it
	# npdf_2 <-  expression(mu*log(x))
	# D2_2 <- D(D(npdf_2,'mu'),'mu')

	D2 <- length(SL) * eval(D2_1,list(a=SLmin,b=SLmax,mu=CI[,1]))

	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)

	int <- qnorm(0.975)*sqrt(1/D2)
	CI[,2] <- CI[,1] - int
	CI[,3] <- CI[,1] + int

  return(CI)
}

#######################################
# As we are assuming that the step length and turning angles are independent
# from one another in the BW, TBW, CRW, TCRW, 
# we can estimate the the CI for the lambdas and kappas independtly.
# Thus I am not having seperate function for BW, TBW, CRW. TCRW
# But rather CI.E, CI.TE, CI.K
# CI.E for lambda
# CI.TE for truncated lambda
# CI.k for kappa


#######################################
# E
# Based on nll.BW

ciE <- function(SL,TA,SLmin,mleM){
  # Table for CI
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimate
  CI[,1] <- mleM[1]
	
  parF <- list('SLmin'=SLmin)

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# In the case of the BW
	# there is an analytical solution for the second derivative
	# D2 = -n/lambda^2
	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	# 1/sqrt(-n/lambda^2) = lambda/sqrt(-n)
	int <- qnorm(0.975)*CI[,1]/sqrt(length(SL))
	CI[,2] <- CI[,1] - int
	CI[,3] <- CI[,1] + int

	return(CI)
}


#######################################
# K
# This is actually based on nll.CRW

ciK <- function(SL,TA,SLmin,mleM,graph=TRUE){
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[2])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  CI[1] <- mleM[2]

  parF <- list('SLmin'=SLmin, 'lambda'=mleM[1])
	
	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.

	
	##
	# Kappa
	# The LL of the vonmises distribution when a=0 is
	# exp(kapp*cos(SL))/(2*pi*besselI(kapp,0))
	# See Forbes et al. 2011
	# nu is the order -> in this case 0
	# I couldn't find the 2nd derviative with R
	# SO I did it with Mathematica
	D2_npdf <- expression((besselI(k,0)^2 - 2* besselI(k,1)^2 + besselI(k,0)*besselI(k,2))/(2*besselI(k,0)^2))
	D2 <- length(TA)*eval(D2_npdf,list(k=CI[,1]))

	int <- qnorm(0.975)*sqrt(1/D2)
	CI[2] <- CI[,1] - int
	CI[3] <- CI[,1] + int

	return(CI)
}

#######################################
# TE
# Based on nll.TBW

ciTE <- function(SL,TA,SLmin,SLmax,mleM){
  # Table for CI
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  
  # Parameter estimate
  CI[1] <- mleM[1]
  
	parF <- list('SLmin'=SLmin,'SLmax'=SLmax)

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.

	# In the case of the TBW
	# there is an analytical solution for the second derivative
	# But it is ugly.
	# So I'm using R to get the second derivative
	# I first deriving the negative log pdf in 2 parts
	# Part 1: part without SL
	npdf_1 <- expression(-(log(lambda) - log(exp(-lambda*a)-exp(-lambda*b))))
	D2_1 <- D(D(npdf_1,'lambda'),'lambda')
	# Part 2: part with SL, but D2_2 = 0
	# SO I exclude
	# npdf_2 <- expression(lambda*x)
	# D2_2 <- D(D(npdf_2,'lambda'),'lambda')

	D2 <- length(SL)*eval(D2_1,list(lambda=CI[,1], a=SLmin, b=SLmax))

	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	int <- qnorm(0.975)*sqrt(1/D2)
	CI[2] <- CI[,1] - int
	CI[3] <- CI[,1] + int

	return(CI)
}
