#######################################
# Created By: Marie Auger-Methe
# Date created: November 21, 2011
# Updated: January 26, 2012

# This script contains the functions for computing the profile likelihood (PL) CI
# and the quadratic approximation CI (using the Hessian).
# The fx for the quad. app. CI CI.Hessian() can be used directly
# to calculate the CI for all parameters estimated through mle of a model
# For the PL CI you need the overall function that calculate the PL CI (CI.PL())
# for a given parameter,
# and the individual fx for each model that uses CI.PL()
# to calculate the PL CIs for all the parameters estimated through mle
# Note that the range of parameter values for which the PL is calculated is based
# on the quad. app. CI.

#######################################
CI.Hessian <- function(SL,TA_N, parMLE, trans.par, M, parF, NLL){

	# Create matrix for Hessian C.I. results
	parCI <- matrix(NA,nrow=length(parMLE),ncol=2)
	colnames(parCI) <- c("L95CI", "U95CI")
	rownames(parCI) <- names(parMLE)

	# Get a numerial estimate of the Hessian from the optimisation routine
	# but this is for the working parameter
	# the nlm often crashes for CCRW_HMM so can't use it for hessian
	H <- hessian(NLL, trans.par(parMLE),SL=SL,TA_N=TA_N,parF=parF)

	# Transform for real parameter
	# From Zucch & MacD 2009
	# G^{-1} = M' %*% H^{-1} %*% M
	# Where H is the hessian of the working paramter
	# G is the hessian of the true parameters
	# M is define by m_ij = d theta_j / d phi_i
	# (where theta is true parameter and phi are working parameter)

	# although I feel like there is an easier solution
	# For X: d plogis(phi) /d phi
	# exp(phi)/(1+exp(phi))^2
	# For lambda_F, lambda_T, kappa_T: (doesn't include truncation constraint?!?)
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


		# To get the corelation between the parameters
		# divide sigma_ab/(sigma_a*sigma_b)
	#	cor_par <- cov2cor(G_inv) 
	
		# According to Ben Bolker's book:
		# N(alpha) * sqrt(1/hessian)
		# So I think N(alpha) * sqrt(diag(G_inv))

		int <- qnorm(0.975)*sqrt(diag(G_inv))
	
		parCI[,1] <- parMLE - int
		parCI[,2] <- parMLE + int
	}
	return(parCI)	

}

#######################################
CI.PL.EM <- function(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mnll.m, notMisLoc, B=100,graph=T){

#######################################
# This function calculates the 95% confidence interval
# using the the profile likelihood method
# described in Bolker's book

# Edwards 2008 and Edwards 2011 use the profile likelihood method for the CI
# According to Bolker's book it is the most robust method to get CI
# Unlike the quadratic approximation or the slice likelihood
# it is less affected by correlation between parameters and
# departure from quadractic shape of the likleihood surface at the mle
# But profile likelihood CI are slow to compute.

# I'm using the mle as the starting value for the minimisation
# rather than using multiple starting values.
# The problem with the likelihood profile is that it doesn't work
# when the value is on the edge of the allowable range
# This is similar to the hessian values
# and I think this is why Zucchini & MacDonald use the parametric bootstrap.
# But parametric boothstrap is really slow

# Function inputs
# SL: step length
# TA: relative turn angle
# parI: parameter of interest (for which the CI is computed with this fx) (not transformed)
# parMLE.t: MLE values for all the parameters that are estimated through MLE (already transformed)
# NLL: neg LL function for the model (one of the fx in nll.fx.R)
# B: number of values of parI for which the min neg LL is estimated
# rang.b: range boundaries for parI
# mnll.m: min neg LL of the model when all parameters are minimised

	# When the mle is at the boundary of the possible values it can take
	# you can run into problems where there is such a big discrepancy between
	# the lower and upper interval range
	# that equaly dividing the range might result in problems
	# to somewhat go around this problem
	# I divide each branch equally
	mnll_i <- floor(B/2)
	rang.parI <- numeric(B)
	rang.parI[1:mnll_i] <- seq(rang.b[1],parI,length=mnll_i)
	rang.parI[mnll_i:B] <- seq(parI,rang.b[2],length=(B-mnll_i+1))
	names(rang.parI) <- rep(names(parI),B)

	# Set memory for loop results
	rang.mnll <- numeric(B)
	mnllRes <- vector("list", length=5)
 
	for (i in 1:B){
		mnllRes <- EM_CCRW_HMM_CI(SL,TA_N,parMLE,missL,notMisLoc,parF=c(parF,rang.parI[i]))
		rang.mnll[i] <- mnllRes$mllk
		mnllRes$message
	}

	# Step 2. Divide into two branch (Lower and upper)
	# Since we divided around the mle we can use the index from above
	prof.lower <- rang.mnll[1:mnll_i]
	prof.l.avec <- rang.parI[1:mnll_i]

	prof.upper <- rang.mnll[mnll_i:B]
	prof.u.avec <- rang.parI[mnll_i:B]

	# Step 3. Using linear interpolation to get the value of X for which
	# neg LL (for the parameter range) =
	# neg LL (from the real MLE) + chisqdist(0.95)/2

	thr_95 <- mnll.m + qchisq(0.95,1)/2

  # So we get to good peak
  if(any(prof.lower > thr_95, na.rm=T)){
    im <- max(which(prof.lower > thr_95))
    if(im >1){
      prof.lower[1:(im-1)] <- NA 
    }
  }
  if(any(prof.upper > thr_95, na.rm=T)){
    # So we get to good peak
    im <- min(which(prof.upper > thr_95))
    if(im < length(prof.upper)){
      prof.upper[(im+1):length(prof.upper)] <- NA 
    }
  }

	lci <- approx(prof.lower,prof.l.avec,
		xout = thr_95)

	uci <- approx(prof.upper,prof.u.avec,
		xout = thr_95, ties="ordered")

	res <- c('LCI'=lci$y, 'UCI'=uci$y)

	if(graph==T){
		# The likelihood profile is pretty ugly for the values at the boundary
		windows()
		plot(rang.mnll~rang.parI,
			main=paste("Likelihood profile of",names(parI)),
			xlab="Parameter value", ylab="Neg LL")
		abline(h=thr_95)
		abline(v=res)
	}

	return(res)
}

#######################################
CI.Slice <- function(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mnll.m, B=100,graph=T){

#######################################
# This function calculates the 95% confidence interval
# using the the profile likelihood method
# described in Bolker's book
# But this is for models that have 1 parameter estimted with the mle
# So really in this case it is either the true likelihood curve
# or a slice of the likelihood if the other paramters are estimated
# with another method
# e.g. SLmin is the min SL

# The problem with the likelihood profile is that it doesn't work
# when the value is on the edge of the allowable range
# I think this is why Zucchini & MacDonald use the parametric bootstrap.
# But parametric boothstrap is really slow

# Function inputs
# SL: step length
# TA: relative turn angle
# parI: parameter of interest (for which the CI is computed with this fx) (not transformed)
# NLL: neg LL function for the model (one of the fx in nll.fx.R)
# B: number of values of parI for which the min neg LL is estimated
# rang.b: range boundaries for parI
# mnll.m: min neg LL of the model when all parameters are minimised

	# When the mle is at the boundary of the possible values it can take
	# you can run into problems where there is such a big discrepancy between
	# the lower and upper interval range
	# that equaly dividing the range might result in problems
	# to somewhat go around this problem
	# I divide each branch equally
	mnll_i <- floor(B/2)
	rang.parI <- numeric(B)
	# In some cases, expecially for mleof k
	# which is bounded at 0
	# the parI is the same as rang.b[1]
	# which causes problem and I have to return to orginal seq

	if(rang.b[1]==parI){
		rang.parI <- seq(rang.b[1],rang.b[2],length=B)
	}else{
		rang.parI[1:mnll_i] <- seq(rang.b[1],parI,length=mnll_i)
		rang.parI[mnll_i:B] <- seq(parI,rang.b[2],length=(B-mnll_i+1))
	}
	names(rang.parI) <- rep(names(parI),B)

	# Creating a matrix that will save the different minimiztion results
	rang.mnll <- numeric(B)

	for (i in 1:B){
		rang.mnll[i] <- NLL(SL=SL,TA_N=TA_N,trans.par(rang.parI[i]),parF=parF)
	}

  if(rang.b[1]==parI){
    # Because the mle value is at the lower boundary you can only get an upper CI
    
    # Lower CI set as NA
    lci <- list(y=NA)
    
    # Get upper CI
    uci <- approx(rang.mnll, rang.parI,
                  xout = mnll.m + qchisq(0.95,1)/2)
    
  }else{
    # Step 2. Divide into two branch (Lower and upper)
    # Since we divided around the mle we can use the index from above
    prof.lower <- rang.mnll[1:mnll_i]
    prof.l.avec <- rang.parI[1:mnll_i]
    
    prof.upper <- rang.mnll[mnll_i:B]
    prof.u.avec <- rang.parI[mnll_i:B]
    
    # Step 3. Using linear interpolation to get the value of X for which
    # neg LL (for the parameter range) =
    # neg LL (from the real MLE) + chisqdist(0.95)/2
    
    lci <- approx(prof.lower,prof.l.avec,
                  xout = mnll.m + qchisq(0.95,1)/2)
    
    uci <- approx(prof.upper,prof.u.avec,
                  xout = mnll.m + qchisq(0.95,1)/2)
  }


	res <- c('LCI'=lci$y, 'UCI'=uci$y)

	if(graph==T){
		# The likelihood profile is pretty ugly for the values at the boundary
		windows()
		plot(rang.mnll~rang.parI,
			main=paste("Likelihood profile of",names(parI)),
			xlab="Parameter value", ylab="Neg LL")
		abline(h=(mnll.m + qchisq(0.95,1)/2))
		abline(v=res)
	}

	return(res)
}

# Warning message used in CI.CCRW_IM & CI.CCRW_HMM
# If you can't calculate the Hessian CI
# You can automatically set the range for the PL
# Need to do it manually
warn.mess <- 'Hessian CI gives 0, NA, NaN, or any non-finite values. 
	Cannot used the Hessian CI to automatically set the range of values
	for which the PL is calculated.
	Need to manually set range for PL CI, using CI.PL().'

#######################################
# CCRW_HMM

CI.CCRW_HMM <- function(SL,TA_N,SLmin,missL,notMisLoc,mleM,graph=T){
	NLL <- nll.CCRW_HMM
	parMLE <- mleM[1:5]
	parF <- list('SLmin'=SLmin)
	trans.par <- function(x){
		x[1] <- qlogis(x[1])
		x[2] <- qlogis(x[2])
		x[3] <- log(x[3] - .Machine$double.xmin)
		x[4] <- log(x[4] - .Machine$double.xmin)
		x[5] <- log(x[5])
		return(x)
	}

	delta_F <- mleM[6]
	delta_T <- mleM[7]

	# According to Zucchini and MacDonald (2009)
	M <- diag(c(
		parMLE[1]*(1-parMLE[1]),
		parMLE[2]*(1-parMLE[2]),
		parMLE[3],parMLE[4],parMLE[5]))

	# Note that although I put places for the deltas
	# I won't get any CIs
	CI <- matrix(NA,nrow=length(parMLE),ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_H", "U95CI_H","L95CI_PL", "U95CI_PL")

	CI[,1] <- parMLE

  if(any(is.na(mleM[1:7]))){
	  warning("The EM-algorithm gives NA values for the parameters")
	}else{
	  CI[,2:3] <- CI.Hessian(SL,TA_N, parMLE, trans.par, M,
	                         parF=c(parF,delta_F,delta_T,list("missL"=missL)), NLL)
	}
	
	parMLE <- mleM[1:6]


	# There is issue when the Hessian CI gives values that are outside the allowed range
	# Setting it to the boudary itself often creates problems
	# since lambda is 1/mean(SL-SLmin) it should not be smaller than if you use the SLmax
	# To be safe I add 1 order of magnitude
	lambda_min <- 0.1/max(SL-SLmin)
	


	if(any(!is.finite(CI[,2:3])|(CI[,2]-CI[,3]==0))){
		warning(warn.mess)
	}else{
		# a11
		parI <- parMLE[1]
		int <- 3* (parI - CI[1,2])
		rang.b <- c(max(parI-int,0),min(parI + int,0.99999))
		CI[1,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
			notMisLoc, B=100,graph)
    if(any(is.na(CI[1,4:5]))){
      if(is.na(CI[1,4])){
        rang.b[1] <- 0
      }
      if(is.na(CI[1,5])){
        rang.b[2] <- 1
      }
      dev.off(dev.cur())
      CI[1,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
                            notMisLoc, B=150,graph)
      if(any(is.na(CI[1,4:5]))){
        if(is.na(CI[1,4])){
          CI[1,4] <- 0
        }
        if(is.na(CI[1,5])){
          CI[1,5] <- 1
        }
      }
    }
		if(graph==T){abline(v=CI[1,2:3],col='red')}

		# a22
		parI <- parMLE[2]
		int <- 3* (parI - CI[2,2])
		rang.b <- c(max(parI-int,0),min(parI + int,0.99999))
		CI[2,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
		                      notMisLoc, B=100,graph)
		if(any(is.na(CI[2,4:5]))){
		  if(is.na(CI[2,4])){
		    rang.b[1] <- 0
		  }
		  if(is.na(CI[2,5])){
		    rang.b[2] <- 1
		  }
		  dev.off(dev.cur())
		  CI[2,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
		                        notMisLoc, B=150,graph)
		  if(any(is.na(CI[2,4:5]))){
		    if(is.na(CI[2,4])){
		      CI[2,4] <- 0
		    }
		    if(is.na(CI[2,5])){
		      CI[2,5] <- 1
		    }
		  }
		}
		if(graph==T){abline(v=CI[2,2:3],col='red')}


		# lambda_F
		parI <- parMLE[3]
		int <- 3* (parI - CI[3,2])
		rang.b <- c(max(parI-int,lambda_min),parI + int)
		CI[3,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
			notMisLoc,B=150,graph)
		if(any(is.na(CI[3,4:5]))){
      int <- 10*int
		  if(is.na(CI[3,4])){
		    rang.b[1] <- max(parI-int,lambda_min)
		  }
		  if(is.na(CI[3,5])){
		    rang.b[2] <- parI + int
		  }
		  dev.off(dev.cur())
		  CI[3,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
		                        notMisLoc, B=150,graph)
		}
		if(graph==T){abline(v=CI[3,2:3],col='red')}
	
		# lambda_T
		parI <- parMLE[4]
		int <- 3* (parI - CI[4,2])
		rang.b <- c(max(parI-int,lambda_min),parI + int)
		CI[4,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
			notMisLoc, B=150,graph)
		if(any(is.na(CI[4,4:5]))){
		  int <- 10*int
		  if(is.na(CI[4,4])){
		    rang.b[1] <- max(parI-int,lambda_min)
		  }
		  if(is.na(CI[4,5])){
		    rang.b[2] <- parI + int
		  }
		  dev.off(dev.cur())
		  CI[4,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'], 
		                        notMisLoc, B=150,graph)
		}
		if(graph==T){abline(v=CI[4,2:3],col='red')}

		# kappa_T
		parI <- parMLE[5]
		int <- 4*(parI - CI[5,2])
    int <- if(int < 1.5){1.5}else{int}
		#int <- 2*(parI - CI[5,2])
		#rang.b <- c(max(parI-int,0),parI + int)
		rang.b <- c(0,parI + int)
		CI[5,4:5] <- CI.PL.EM(SL,TA_N, parI, parMLE, missL, parF, NLL, rang.b, mleM['mnll'],
			notMisLoc, B=100,graph)
		if(is.na(CI[5,4])){
		  CI[5,4] <- 0
		}
		if(graph==T){abline(v=CI[5,2:3],col='red')}
	}
	return(CI)
}


#######################################
# LW

CI.LW <- function(SL,TA_N,SLmin,mleM,graph=T){
	NLL <- nll.LW
	parMLE <- mleM[1]
	parF <- list('SLmin'=SLmin)
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# It's particularly useful to find good intervals
	# to do the PL

	# The PL likelihood is more a slice in this case
	# since the other parameter estimated: a
	# is fixed to its value
	CI <- matrix(NA,nrow=length(parMLE),ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_N", "U95CI_N","L95CI_PL", "U95CI_PL")

	CI[,1] <- parMLE

	# In the case of the LW
	# there is an analytical solution for the second derivative
	# D2 = -n/(mu-1)^2
	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	# 1/sqrt(-n/(mu-1)^2) = (mu-1)/sqrt(-n)
	int <- qnorm(0.975)*(parMLE-1)/sqrt(length(SL))
	CI[,2] <- parMLE - int
	CI[,3] <- parMLE + int

	# mu
	parI <- parMLE[1]
	# I use the normal approximation to set the range for the likelihood curve
	int <- 2 * int
	rang.b <- c(max(parI-int,1.000001), parI + int)
	CI[,4:5] <- CI.Slice(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mleM['mnll'], B=100,graph)
	if(graph==T){abline(v=CI[,2:3],col='red')}
	return(CI)
}

#######################################
# TLW

CI.TLW <- function(SL,TA_N,SLmin,SLmax,mleM,graph=T){
	NLL <- nll.TLW
	parMLE <- mleM[1]
	parF <- list('SLmin'=SLmin,'SLmax'=SLmax)
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# It's particularly useful to find good intervals
	# to do the PL

	# The PL likelihood is more a slice in this case
	# since the other parameters estimated: SLmin & SLmax
	# are fixed to their estimated values

	CI <- matrix(NA,nrow=length(parMLE),ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_N", "U95CI_N","L95CI_PL", "U95CI_PL")

	CI[,1] <- parMLE

	# In the case of the TLW
	# there is an analytical solution for the second derivative
	# But it is ugly.
	# So I'm using R to get the second derivative
	# I first deriving the negative log pdf in two parts
	# Part 1 doesn't include the SL 
	npdf_1 <- expression(-(log(mu-1) - log(a^(1-mu) - b^(1-mu))))
	D2_1 <- D(D(npdf_1,'mu'),'mu')
	# Part 2 includes the SL
	# but npdf_2 = 0 so I am not running it
	# npdf_2 <-  expression(mu*log(x))
	# D2_2 <- D(D(npdf_2,'mu'),'mu')

	D2 <- length(SL) * eval(D2_1,list(a=SLmin,b=SLmax,mu=parMLE))

	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)

	int <- qnorm(0.975)*sqrt(1/D2)
	CI[,2] <- parMLE - int
	CI[,3] <- parMLE + int

	# mu
	parI <- parMLE[1]
	int <- 2 * int
	rang.b <- c(max(parI-int,1.000001), parI + int)
	CI[,4:5] <- CI.Slice(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mleM['mnll'], B=100,graph)
	if(graph==T){abline(v=CI[,2:3],col='red')}
	return(CI)
}

#######################################
# I part to be efficient
# I am estimating the the CI for the lambdas and kappas independtly
# and thus I am not having seperate function for BW, TBW, CRW. TCRW
# But rather CI.E, CI.TE, CI.K
# This is possible because in the models CRW & TCRW,
# where both kappa & lambda needs to be estimated,
# the 2 paremeters to estimate with mle are completely indepedent.
# So the CI can be estimated for each of the parameter seperately.


#######################################
# E
# This is actually based on nll.BW

CI.E <- function(SL,TA_N,SLmin,mleM,graph=T){
	NLL <- nll.BW
	parMLE <- mleM[1]
	parF <- list('SLmin'=SLmin)
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# It's particularly useful to find good intervals
	# to do the PL

	# The PL likelihood is more a slice in this case
	# since the other parameter estimated: SLmin
	# is fixed to its estimated value
	CI <- matrix(NA,nrow=length(parMLE),ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_N", "U95CI_N","L95CI_PL", "U95CI_PL")

	CI[,1] <- parMLE

	# In the case of the BW
	# there is an analytical solution for the second derivative
	# D2 = -n/lambda^2
	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	# 1/sqrt(-n/lambda^2) = lambda/sqrt(-n)
	int <- qnorm(0.975)*parMLE/sqrt(length(SL))
	CI[,2] <- parMLE - int
	CI[,3] <- parMLE + int

	# lambda
	parI <- parMLE[1]

	# I use the normal approximation to set the range for the likelihood curve
	int <- 2 * int
	rang.b <- c(max(parI-int,1e-100),parI + int)
	CI[,4:5] <- CI.Slice(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mleM['mnll'], B=100,graph)
	if(graph==T){abline(v=CI[,2:3],col='red')}
	return(CI)
}


#######################################
# K
# This is actually based on nll.CRW

CI.K <- function(SL,TA_N,SLmin,mleM,graph=T){
	NLL <- nll.CRW
	parMLE <- mleM[2]
	parF <- list('SLmin'=SLmin, 'lambda'=mleM[1])
	# Some models used transformed parameters in nll
	# Thus trans.par called in an input of CI.slice
	trans.par <- function(x){x}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# It's particularly useful to find good intervals
	# to do the PL


	# The PL likelihood is more a slice in this case
	# since the other parameter estimated: SLmin
	# is fixed to its estimated value
	CI <- matrix(NA,ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_N", "U95CI_N","L95CI_PL", "U95CI_PL")

	CI[1] <- parMLE
	
	##
	# Kappa
	# The LL of the vonmises distribution when a=0 is
	# exp(kapp*cos(SL))/(2*pi*besselI(kapp,0))
	# See Forbes et al. 2011
	# nu is the order -> in this case 0
	# I couldn't find the 2nd derviative with R
	# SO I did it with Mathematica
	D2_npdf <- expression((besselI(k,0)^2 - 2* besselI(k,1)^2 + besselI(k,0)*besselI(k,2))/(2*besselI(k,0)^2))
	D2 <- length(TA_N)*eval(D2_npdf,list(k=parMLE))

	int <- qnorm(0.975)*sqrt(1/D2)
	CI[2] <- parMLE - int
	CI[3] <- parMLE + int

	parI <- parMLE
	names(parI) <- 'kapp'
	# I use the normal approximation to set the range for the likelihood curve
	#int <- 2 * int
	int <- 3 * int
	int <- if(int < 1){1}else{int}
  #rang.b <- c(max(parI-int,0),parI + int)
	rang.b <- c(0,parI + int)
	CI[4:5] <- CI.Slice(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mleM['mnll'], B=100,graph)
	if(graph==T){abline(v=CI[2:3],col='red')}
	if(is.na(CI[,4])){
	  CI[4] <- 0
	}
	return(CI)
}

#######################################
# TE
# This is actually based on nll.TBW

CI.TE <- function(SL,TA_N,SLmin,SLmax,mleM,graph=T){
	NLL <- nll.TBW
	parMLE <- mleM[1]
	parF <- list('SLmin'=SLmin,'SLmax'=SLmax)
	# TBW uses a transformation in nll.TWB
	trans.par <- function(x){log(x)}

	# Similar to using the hessian
	# we can use the second derivate of the likelihood
	# to estimate the normal approximation of the C.I.
	# It's particularly useful to find good intervals
	# to do the PL

	# The PL likelihood is more a slice in this case
	# since the other parameters estimated: SLmin & SLmax
	# are fixed to their estimated values
	CI <- matrix(NA,ncol=5)
	rownames(CI) <- names(parMLE)
	colnames(CI) <- c("estimate","L95CI_N", "U95CI_N","L95CI_PL", "U95CI_PL")

	CI[1] <- parMLE

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

	D2 <- length(SL)*eval(D2_1,list(lambda=parMLE,a=SLmin,b=SLmax))

	# The normal approximation is
	# N(1-alpha) * sqrt(1/D2)
	int <- qnorm(0.975)*sqrt(1/D2)
	CI[2] <- parMLE - int
	CI[3] <- parMLE + int

	# lambda
	parI <- parMLE

	# Unlike the other models which looks for the mle
	# of only one parameter
	# The TBW a uses transformed parameter
	# so I added trans.par

	# I use the normal approximation to set the range for the likelihood curve
	int <- 2 * int
	rang.b <- c(max(parI-int,1e-100),parI + int)
	CI[4:5] <- CI.Slice(SL,TA_N, parI, parF, trans.par, NLL, rang.b, mleM['mnll'], B=100,graph)
	if(graph==T){abline(v=CI[2:3],col='red')}
	return(CI)
}