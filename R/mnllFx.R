#######################################
# Maximum likelihood estimates
# Finding the minimum neg LL

# When I'm optimizing multiple parameters
# I'm using multiple sets of starting values
# to be sure to get the global minimum and a local minimum.
# I'm also using the data to find appropriate starting values.

# Because we are assuming that the TA is independent of the SL
# we can minimise the neg LL of each element seperately
# for the BW, CRW, TBW, TCRW, LW, TLW.
# Not for the CCRWs because the probability of being in each
# behaviour depends on both the TA and the SL.

################################
# CCRW

mnllCCRW <- function(SL, TA, TA_C, missL, notMisLoc, SLmin, tol=5e-5){

	######################################################################
	# Setting starting values
	# lI & lE
  # row 1:  starting value of lI based on the 1st (shortest) quantile of SL-SLmin
  #         and lE on the 3rd (longuest) quantile of SL-SLmin
  # row 2:  starting value of lI and lE are based on the median SL-SLmin
	l0 <- matrix(1/quantile((SL-SLmin), c(0.25, 0.50, 0.75, 0.50)),ncol=2)
	nL <- length(l0[,1])

	# According to Zucchini and MacDonald 2009,
	# starting values for the transition probability matrix should be symmetric (gII=gEE).
	# So for highest value we use symetric probability.
  # Although the probability of remaining in the intensive behavior should always be
  # high those of the extensive behavior can be low. So we explore low probability values.
	g0 <- c(0.1,0.3,0.5,0.7,0.9)
	nA <- length(g0)

	# Because we use symetric transition probabilty matrices, 
	# the stationary distribution of the behaviours are 0.5.
	# Therefore k0 can be based on the 0.5 quantile
	TA_T <- TA_C[TA>quantile(TA, 0.25) & TA<quantile(TA, 0.75)]
	k0 <- mle.vonmises(TA_T, mu=circular(0))$kappa
  
	par0 <- cbind(rep(g0,each=nL),rep(g0,each=nL),rep(l0[,1],nA),rep(l0[,2],nA),k0)

	# Creating a matrix that will save the minimiztion results
	mnll <- matrix(NA, ncol=9, nrow=nrow(par0))
	colnames(mnll) <- c("gII", "gEE", "lI", "lE", "kE", "dI", "dE", "a", "mnll")
  mnll[,'a'] <- SLmin

	D_0 <- c(0.5,0.5) # First true states
	
	gamm <- function(x){
		matrix(c(x[1],(1-x[1]),(1-x[2]),x[2]), byrow=TRUE, nrow=2)
	}

	for (i in 1:nrow(par0)){
		mnllRes <- tryCatch(emHMM(SL,TA,missL,SLmin,
			par0[i,3:4],gamm(par0[i,1:2]),D_0,par0[i,5],notMisLoc,tol=tol),
			error=function(e) list('gamma'=NA,'lambda'=NA,
			'kapp'=NA,'delta'=NA,'mllk'=NA))
		mnll[i,1:2] <- mnllRes$gamma[c(1,4)]
		mnll[i,3:4] <- mnllRes$lambda
		mnll[i,5] <- mnllRes$kapp
		mnll[i,6:7] <- mnllRes$delta
		mnll[i,'mnll'] <- mnllRes$mllk

	}

	mnll <- mnll[which.min(mnll[,'mnll']),]

	if (length(mnll)==0){ # In case no minimization was able to get good values
		mleCCRW <-rep(NA,13)
		names(mleCCRW) <- c("gII", "gEE", "lI", "lE", "kE",
			"dI", "dE", "a", "mnll", "I*","E*", "AIC", "AICc")
	}else{

		# For the comparison purposes and for fitting to observed data
		# I'm estimating the stationary distribution of the HMM
		e_R <- eigen(t(gamm(mnll[1:2]))) # Get the eigenvalues and eigenvectors
		delta <- e_R$vectors[,1] # To get the dominant eigen vector
		delta <- delta/(sum(delta)) # Normalised dominant eigen vector
		names(delta) <- c("I*","E*")

    # According to Burnham and Anderson (2002)
    # AIC = 2*nll +2*k
    # AICc = AIC + 2*k*(k+1)/(n-K-1)
		# Number of parameters estimated, k=7:
    #     2 lambdas + 1 delta + 2 trans prob + 1 kappa + SLmin
    AICCCRW <- matrix(NA, ncol=2)
		names(AICCCRW) <- c("AIC", "AICc")
		AICCCRW[,1] <- 14 + 2*mnll["mnll"]
		AICCCRW[,2] <- AICCCRW[,1] + 112/(length(SL)-8)
		mleCCRW <- c(mnll, delta, AICCCRW)
	}
	
	return(mleCCRW)
}


#########
# Using direct numerical minimization (instead of EM-algorithm)
mnllCCRWdn <- function(SL, TA, TA_C, missL, SLmin){
  
  ######################################################################
  # Setting starting values
  # lI & lE
  # row 1:  starting value of lI based on the 1st (shortest) quantile of SL-SLmin
  #         and lE on the 3rd (longuest) quantile of SL-SLmin
  # row 2:  starting value of lI and lE are based on the median SL-SLmin
  l0 <- matrix(1/quantile((SL-SLmin), c(0.25, 0.50, 0.75, 0.50)),ncol=2)
  nL <- length(l0[,1])
  
  # According to Zucchini and MacDonald 2009,
  # starting values for the transition probability matrix should be symmetric (gII=gEE).
  # So for highest value we use symetric probability.
  # Although the probability of remaining in the intensive behavior should always be
  # high those of the extensive behavior can be low. So we explore low probability values.
  g0 <- c(0.1,0.3,0.5,0.7,0.9)
  nA <- length(g0)
  
  # Because we use symetric transition probabilty matrices, 
  # the stationary distribution of the behaviours are 0.5.
  # Therefore k0 can be based on the 0.5 quantile
  TA_T <- TA_C[TA>quantile(TA, 0.25) & TA<quantile(TA, 0.75)]
  k0 <- mle.vonmises(TA_T, mu=circular(0))$kappa
  
  par0 <- cbind(qlogis(rep(g0,each=nL)),qlogis(rep(g0,each=nL)),log(rep(l0[,1],nA)),log(rep(l0[,2],nA)),log(k0))
  
  # Creating a matrix that will save the minimiztion results
  mnll <- matrix(NA, ncol=7, nrow=nrow(par0))
  colnames(mnll) <- c("gII", "gEE", "lI", "lE", "kE", "a", "mnll")
  mnll[,'a'] <- SLmin
  
  for (i in 1:nrow(par0)){
    mnllRes <- tryCatch(optim(par0[i,],nllCCRWdn,SL=SL,TA=TA, parF=list(SLmin=min(SL),missL=missL)),
                        error=function(e) list(par=rep(NA,5),'value'=NA))
    mnll[i,1:2] <- plogis(mnllRes$par[1:2])
    mnll[i,3:4] <- .Machine$double.xmin + exp(mnllRes$par[3:4])
    mnll[i,5] <- exp(mnllRes$par[5])
    mnll[i,'mnll'] <- mnllRes$value
  }
  
  mnll <- mnll[which.min(mnll[,'mnll']),]
  gamm <- function(x){
    matrix(c(x[1],(1-x[1]),(1-x[2]),x[2]), byrow=TRUE, nrow=2)
  }
  
  if (length(mnll)==0){ # In case no minimization was able to get good values
    mleCCRW <- rep(NA,11)
    names(mleCCRW) <- c("gII", "gEE", "lI", "lE", "kE",
                        "a", "mnll", "I*","E*", "AIC", "AICc")
  }else{
    
    # For the comparison purposes and for fitting to observed data
    # I'm estimating the stationary distribution of the HMM
    e_R <- eigen(t(gamm(mnll[1:2]))) # Get the eigenvalues and eigenvectors
    delta <- e_R$vectors[,1] # To get the dominant eigen vector
    delta <- delta/(sum(delta)) # Normalised dominant eigen vector
    names(delta) <- c("I*","E*")
    
    # According to Burnham and Anderson (2002)
    # AIC = 2*nll +2*k
    # AICc = AIC + 2*k*(k+1)/(n-K-1)
    # Number of parameters estimated (1 less because no delta, using stationary distribution), k=6:
    #     2 lambdas + 2 trans prob + 1 kappa + SLmin
    AICCCRW <- matrix(NA, ncol=2)
    names(AICCCRW) <- c("AIC", "AICc")
    AICCCRW[,1] <- 12 + 2*mnll["mnll"]
    AICCCRW[,2] <- AICCCRW[,1] + 84/(length(SL)-7)
    mleCCRW <- c(mnll, delta, AICCCRW)
  }
  return(mleCCRW)
}

#########
# Using direct numerical minimization (instead of EM-algorithm)
mnllCCRWww <- function(SL, TA, TA_C, missL){
  
  ######################################################################
  # Setting starting values
  # scI & scE
  # based on the smallest quantiles
  sc0 <- matrix(quantile(SL, c(0.15, 0.25, 0.85, 0.75)),ncol=2)
  nL <- length(sc0[,1])
  # Shape paremeter equivalent to exponentional
  sh0 <- rep(1,2)
  
  # According to Zucchini and MacDonald 2009,
  # starting values for the transition probability matrix should be symmetric (gII=gEE).
  # So for highest value we use symetric probability.
  # Although the probability of remaining in the intensive behavior should always be
  # high those of the extensive behavior can be low. So we explore low probability values.
  g0 <- c(0.1,0.3,0.5,0.7,0.9)
  nA <- length(g0)
  
  # Because we use symetric transition probabilty matrices, 
  # the stationary distribution of the behaviours are 0.5.
  # Therefore k0 can be based on the 0.5 quantile
  TA_T <- TA_C[TA>quantile(TA, 0.25) & TA<quantile(TA, 0.75)]
  r0 <- mle.wrappedcauchy(TA_T, mu=circular(0))$rho
  
  par0 <- cbind(qlogis(rep(g0,each=nL)),qlogis(rep(g0,each=nL)),
                log(rep(sc0[,1],nA)),log(rep(sc0[,2],nA)),
                log(sh0[1]),log(sh0[2]),qlogis(r0))
  
  # Creating a matrix that will save the minimiztion results
  mnll <- matrix(NA, ncol=8, nrow=nrow(par0))
  colnames(mnll) <- c("gII", "gEE", "scI", "scE", "shI", "shE", "rE", "mnll")
  
  for (i in 1:nrow(par0)){
    mnllRes <- tryCatch(optim(par0[i,],nllCCRWww,SL=SL,TA=TA, parF=list(missL=missL)),
                        error=function(e) list(par=rep(NA,7),'value'=NA))
    mnll[i,1:2] <- plogis(mnllRes$par[1:2])
    mnll[i,3:6] <- .Machine$double.xmin + exp(mnllRes$par[3:6])
    mnll[i,7] <- plogis(mnllRes$par[7])
    mnll[i,'mnll'] <- mnllRes$value
  }
  mnll <- mnll[which.min(mnll[,'mnll']),]
  
  gamm <- function(x){
    matrix(c(x[1],(1-x[1]),(1-x[2]),x[2]), byrow=TRUE, nrow=2)
  }
  
  if (length(mnll)==0){ # In case no minimization was able to get good values
    mleCCRW <- rep(NA,12)
    names(mleCCRW) <- c("gII", "gEE", "scI", "scE", "shI", "shE", "rE", "mnll", "I*","E*", "AIC", "AICc")
  }else{
    
    # For the comparison purposes and for fitting to observed data
    # I'm estimating the stationary distribution of the HMM
    stdist <- function(gg){
      e_R <- eigen(t(gamm(gg))) # Get the eigenvalues and eigenvectors
      delta <- e_R$vectors[,1] # To get the dominant eigen vector
      delta <- delta/(sum(delta)) # Normalised dominant eigen vector
      names(delta) <- c("I*","E*")
      return(delta)
    }
    delta <- stdist(mnll[1:2]) # Get the eigenvalues and eigenvectors
    
    # According to Burnham and Anderson (2002)
    # AIC = 2*nll +2*k
    # AICc = AIC + 2*k*(k+1)/(n-K-1)
    # Number of parameters estimated (1 less because no delta, using stationary distribution and no SLmin),
    # but extra shape parameter for each weibull step length dist: k=7:
    #     2 scale + 2 shape + 2 trans prob + 1 rho
    AICCCRW <- matrix(NA, ncol=2)
    names(AICCCRW) <- c("AIC", "AICc")
    AICCCRW[,1] <- 14 + 2*mnll["mnll"]
    AICCCRW[,2] <- AICCCRW[,1] + 112/(length(SL)-8)
    mleCCRW <- c(mnll, delta, AICCCRW)
  }
  return(mleCCRW)
}

################################
# LW

mnllLW <- function(SL, TA, SLmin=min(SL)){

	mu <- 1 - 1/(log(SLmin)-mean(log(SL)))
	mnll <- nllLW(SL,TA,mu,parF=list('SLmin'=SLmin))

	# Constrained mu between 1 and 3 if mu_LW is analytical mle is outside
	if(mu <= 1 | mu >3){
		mnll <- optimize(nllLW, interval=c(1,3),SL=SL,TA=TA,parF=list('SLmin'=SLmin))
		mu <- mnll$minimum
		mnll <- mnll$objective
	}
	
	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Number of parameters k=2:
  #     mu + SLmin
	AICLW <- 4 + 2*mnll
	AICcLW <- AICLW + 12/(length(SL)-3)
  
	mleLW <- c(mu, SLmin, mnll, AICLW, AICcLW)
	names(mleLW) <- c("mu", "a", "mnll", "AIC", "AICc")

	return(mleLW)
}

################################
# TLW

mnllTLW <- function(SL,TA,SLmin=min(SL),SLmax=max(SL), conts=TRUE){
	# The data does not give you information on the upper bound of 
	# the LW other than the minimum value this upper bound is 
	# which correspound to the largest step length.
	# Note that there is no analytical solution for the mle of TLW
	# and that mu can =1 see Edwards 2011 supp. (But I need to change the likelihood for it)

  # If contsraining the values to be between 1 & 3, 
  # which are the values consistent with the levy walk hypothesis.
	if(conts){ 
	  mnll <- optimize(nllTLW, interval=c(1,3),SL=SL,TA=TA,parF=list('SLmin'=SLmin,'SLmax'=SLmax))  
	}else{
	  mnll <- optimize(nllTLW, interval=c(-50,50), SL=SL,TA=TA,parF=list('SLmin'=SLmin,'SLmax'=SLmax)) 
	}
 
	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Numbers of parameters k=3:
  #     mu + SLmin + SLmax
	AICTLW <- 6 + 2*mnll[[2]]
  AICcTLW <- AICTLW + 24/(length(SL)-4)
  
	mleTLW <- c(mnll[[1]], SLmin, SLmax, mnll[[2]], AICTLW, AICcTLW)
	names(mleTLW) <- c("mu", "a", "b","mnll", "AIC", "AICc")
	
	return(mleTLW)
}

################################
# BW

mnllBW <- function(SL,TA,SLmin=min(SL)){
	# The only parameter to estimate, lambda, has an analytical solution.
	lambda <- 1/mean(SL-SLmin)
	mnll <- nllBW(SL,TA,lambda,parF=list('SLmin'=SLmin))

	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Numbers of parameter k=2:
  #     lambda + SLmin
	AICBW <- 4 + 2*mnll
  AICcBW <- AICBW + 12/(length(SL)-3)

	mleBW <-  c(lambda,SLmin,mnll,AICBW, AICcBW)
	names(mleBW) <- c("lambda","a","mnll","AIC", "AICc")
	return(mleBW)
}

################################
# CRW

mnllCRW <- function(SL,TA_C, TA,SLmin=min(SL),lambda=(1/mean(SL-SLmin))){
	# Lambda is the same as lambda_BW

	kapp <- mle.vonmises(TA_C,mu=circular(0))$kappa
	mnll <- nllCRW(SL,TA,kapp,lambda,parF=list('SLmin'=SLmin))

	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Numbers of parameters k=3:
  #     lambda + kappa + SL_min
	AICCRW <- 6 +2*mnll
	AICcCRW <- AICCRW + 24/(length(SL)-4)

  mleCRW <- c(lambda,kapp, SLmin, mnll, AICCRW, AICcCRW)
	names(mleCRW) <- c("lambda", "kappa", "a","mnll","AIC", "AICc")
	return(mleCRW)
}

################################
# TBW

mnllTBW <- function(SL,TA,SLmin=min(SL),SLmax=max(SL),lambda_BW=(1/mean(SL-SLmin))){
	# I'm choosing to use nlm(), because optimize require bounds for the parameter.
	# The upper bound of lambda is difficult to asses.
	# I'm using lambda_BW as the starting value.

	mnll <- nlm(nllTBW, log(lambda_BW),SL=SL,TA=TA,parF=list('SLmin'=SLmin,'SLmax'=SLmax))

	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Numbers of parameters k=3:
	#     lambda + SLmin + SLmax
	AICTBW <- 6 + 2*mnll$minimum
	AICcTBW <- AICTBW + 24/(length(SL)-4)
  
	mleTBW <- c(exp(mnll$estimate), SLmin, SLmax,mnll$minimum, AICTBW, AICcTBW)
	names(mleTBW) <- c("lambda", "a", "b", "mnll", "AIC", "AICc")
	return(mleTBW)	
}

################################
# TCRW

mnllTCRW <- function(SL,TA_C,TA,SLmin=min(SL),SLmax=max(SL),
	lambda=exp(nlm(nllTBW, log(1/mean(SL-SLmin)),SL=SL,TA=TA,SLmin=SLmin,SLmax=SLmax)$estimate),
	kapp=mle.vonmises(TA_C,mu=circular(0))$kappa){
	# Same lambda as TBW
	# same kappa as CRW

	mnll <- nllTCRW(SL,TA,lambda,kapp,SLmin,SLmax)

	# According to Burnham and Anderson (2002)
	# AIC = 2*nll +2*k
	# AICc = AIC + 2*k*(k+1)/(n-k-1)
	# Numbers of parameters k=4: 
  #     lambda_TCRW + kappa_TCRW + SLmin + SLmax
	AICCRW <- 8 +2*mnll
  AICcCRW <- AICCRW + 40/(length(SL)-5)
    
	mleTCRW <- c(lambda, kapp, SLmin, SLmax, mnll, AICCRW, AICcCRW)
	names(mleTCRW) <- c("lambda", "kappa", "a", "b","mnll","AIC", "AICc")
	return(mleTCRW)
}
