#######################################
# This set of functions are defining the 
# negative log likelihood (neg LL) for each model
# They are written in great for facilitating their minimisation
# which is done through the fx within mnll.Fx
# As a general rule we use neg LL, instead of the (pos) LL
# because nlm() and optimize() minize by default
# and AIC also use neg. LL.
# As a general rule,
# I'm calculating the LL of each step
# and summing these to get to overall LL
# Note that for the CCRW_HMM
# I'm not directly minimizing the neg LL
# and I am using an EM algorithm
# The nll.CCRW_HMM here is only use for the C.I

#######################################
# CCRW_HMM
# Note that unlike the other nllFx
# this function is not used for when minimizing the neg LL
# This is only used for to calculate the CI.Hessian

nllCCRW <- function(SL,TA,x,parF){
	# ParF should have SLmin and dI & dE and missL

	# Based on ch3 of Zucchini & MacDonald 2009
	# This is likelihood function that can be numerically maximize with respect to the parameters

	####
	# Basic theory
	# So for a vector of observations x of length n (note that in text n is T)

	# We use the forward probability alpha's
	# alpha_1 = delta %*% P(x_1)
	# delta is the initial distribution of the Markov chain (states) (and should sum to 1)
	# P(x_t) is the diagonal matrix of the state-dependent probability
	# of the observation at time t
	# P(x_t) = diag(p_1(x), ..., p_m(x))

	# alpha_t = alpha_{t-1} %*% Gamma %*% P(x_t) 
	# for t = {2,..,n}
	# Gamma is the transition probability matrix for which the rows add to 1

	# L_T = alpha_n %*% t(1)
	# SO the likelihood is recursive

	###
	# Scaling 
	# Instead of using the log(p+q) = lp + log(1 + exp(lq-lp))
	# where p > q and lp = log(p) and lq = log(q)
	# as we did for the CCRW_IM
	# We use the scaled likelihood
	# based on scaling the vector of forward probability  alpha_t
	# See section 3.2 p.46 of Z&M (2009)
	# The scaling (like the method we used before) reduced the potential
	# under- and over- flow (i.e. getting 0 and Inf)

	# We define the scaling weight w_t such that:
	# phi_t = alpha_t / w_t
	# w_t = alpha_t %*% t(1)

	# I think the idea is that the alpha_t becomes smaller and smaller
	# So the sum of the elements of alpha_t, which is w_t is also smaller
	# This makes phi_t bigger

	# To initialise
	# w_1 = alpha_1 %*% t(1) = delta %*% P(x_1) %*% t(1)
	# phi_1 = alpha_1 / w_1

	# Since:
	# alpha_t = alpha_{t-1} %*% Gamma %*% P(x_t)
	# and we can write simplfy the notation by B_t = Gamma %*% P(x_t)
	# w_t %*% phi_t = w_{t-1} %*% phi_{t-1} %*% B_t
	# So L_T = alpha_T %*% t(1) = w_T %*% phi_T %*% t(1)
	# since phi_T = alpha_T / w_T
	# and w_T = alpha_T %*% t(1)
	# phi_T %*% t(1) = alpha_T %*% t(1) / alpha_T %*% t(1) = 1 (scalar)
	# So L_T = w_T %*% phi_T %*% t(1) = w_T (scalar)

	# since w_t %*% phi_t = w_{t-1} %*% phi_{t-1} %*% B_t
	# and phi_t %*% t(1) = 1 (scalar) (and w_t is a scalar)
	# w_t = w_{t-1} %*% phi_{t-1} %*% B_t  %*% t(1)
	# so 
	# L_T = prod_{t=1}^T (phi_{t-1} %*% B_t  %*% t(1))
	# log L_T = sum_{t=1}^T  log(phi_{t-1} %*% B_t  %*% t(1))

	####
	# Set parameters
	n <- length(SL) # Length of time serie, refered as T in text
	t_1 <- matrix(1,nrow=2) # this is 1' or t(1)

	# Set parameters of state-dependent probabilities of observations
	# Uniform angle distribution for F: mu =0, kappa=0
	# Forward persistance for T: mu=0

	#####
	# Parameters to estimate

	# Note that to be able to use unconstrained optimizers
	# I have transformed the parameters
	# plogis constraint values betwen 0 & 1
	# exp() > 0 (actually, because of rounding >=0)
	gII <- plogis(x[1])
	gEE <- plogis(x[2])
	lI <- .Machine$double.xmin + exp(x[3])
	lE <- .Machine$double.xmin + exp(x[4])
	kE <- exp(x[5]) # Note that kappa can be 0 (uniform)


	##
	# To allow to fix some of the parameters
	# This is need for the profile likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}

	# Gamma = the t.p.m
	# Only estimating gII, gEE,
	# a12 & a21 are (1-gII) and (1-gEE), respectively
	# This assures that the row of the t.p.m sum to 1
	# Which mean that you will end up in one of the 2 states at the next time step
	# with a probability of 1
	Gamma = matrix(c(gII,(1-gII),(1-gEE),gEE), nrow=2, byrow=T)

	# To be consistent with the EM algorithm
	# I am not using the stationary distribution of the Markov Chain
	# The initial distribution of the Markov Chain is estimated
	# Constrained to be between 0 and 1
	delta <- c(dI, dE)

	# P(x_t) is the diag of p_i(x_t)
	P <- function(SL_t,TA_t){
		# State 1: F
		# dexp(SL|lI) * dvonmises(TA|muI, kappa_F)
		pI <-   dexp(SL_t-SLmin,lI) * dvm(TA_t, 0, 0)
		# State 2: T
		# dexp(SL|lE) * dvonmises(TA|muE, kE)
		pE <-   dexp(SL_t-SLmin,lE) * dvm(TA_t, 0, kE)

		result <- diag(c(pI,pE))
	}


	########################
	# Creating the variables that will save values during the recursive
	# calculation of the likelihood
	phi <- matrix(NA, ncol=2, nrow = n)

	#####################
	# Initialising the likelihood

	# The weight and scaled forward probabilities of the initial state
	w_1 <- delta %*% P(SL[1],TA[1]) %*% t_1
	phi[1,] <- (1/w_1) %*% delta %*% P(SL[1],TA[1])
	# the log likelihood of the 1st observation
	LL <- log(w_1)
	for (i in 2:n){
		v <- phi[(i-1),] %*% (Gamma %^% missL[i-1]) %*% P(SL[i],TA[i])
		u <- v %*% t_1
		LL <- LL + log(u)
		phi[i,] <- (1 / u) %*% v 
	}
	return(-LL) # Return loglikelihood
}

#######################################
# CCRW for numerical maximization (instead of EM-algorithm)

nllCCRWdn <- function(SL,TA,x,parF=list(SLmin=min(SL),missL=missL)){
  # ParF should have SLmin and missL
  
  # Based on ch3 of Zucchini & MacDonald 2009
  # This is likelihood function that can be numerically maximize with respect to the parameters
  
  ####
  # Basic theory
  # So for a vector of observations x of length n (note that in text n is T)
  
  # We use the forward probability alpha's
  # alpha_1 = delta %*% P(x_1)
  # delta is the initial distribution of the Markov chain (states) (and should sum to 1)
  # P(x_t) is the diagonal matrix of the state-dependent probability
  # of the observation at time t
  # P(x_t) = diag(p_1(x), ..., p_m(x))
  
  # alpha_t = alpha_{t-1} %*% Gamma %*% P(x_t) 
  # for t = {2,..,n}
  # Gamma is the transition probability matrix for which the rows add to 1
  
  # L_T = alpha_n %*% t(1)
  # SO the likelihood is recursive
  
  ###
  # Scaling 
  # Instead of using the log(p+q) = lp + log(1 + exp(lq-lp))
  # where p > q and lp = log(p) and lq = log(q)
  # as we did for the CCRW_IM
  # We use the scaled likelihood
  # based on scaling the vector of forward probability  alpha_t
  # See section 3.2 p.46 of Z&M (2009)
  # The scaling (like the method we used before) reduced the potential
  # under- and over- flow (i.e. getting 0 and Inf)
  
  # We define the scaling weight w_t such that:
  # phi_t = alpha_t / w_t
  # w_t = alpha_t %*% t(1)
  
  # I think the idea is that the alpha_t becomes smaller and smaller
  # So the sum of the elements of alpha_t, which is w_t is also smaller
  # This makes phi_t bigger
  
  # To initialise
  # w_1 = alpha_1 %*% t(1) = delta %*% P(x_1) %*% t(1)
  # phi_1 = alpha_1 / w_1
  
  # Since:
  # alpha_t = alpha_{t-1} %*% Gamma %*% P(x_t)
  # and we can write simplfy the notation by B_t = Gamma %*% P(x_t)
  # w_t %*% phi_t = w_{t-1} %*% phi_{t-1} %*% B_t
  # So L_T = alpha_T %*% t(1) = w_T %*% phi_T %*% t(1)
  # since phi_T = alpha_T / w_T
  # and w_T = alpha_T %*% t(1)
  # phi_T %*% t(1) = alpha_T %*% t(1) / alpha_T %*% t(1) = 1 (scalar)
  # So L_T = w_T %*% phi_T %*% t(1) = w_T (scalar)
  
  # since w_t %*% phi_t = w_{t-1} %*% phi_{t-1} %*% B_t
  # and phi_t %*% t(1) = 1 (scalar) (and w_t is a scalar)
  # w_t = w_{t-1} %*% phi_{t-1} %*% B_t  %*% t(1)
  # so 
  # L_T = prod_{t=1}^T (phi_{t-1} %*% B_t  %*% t(1))
  # log L_T = sum_{t=1}^T  log(phi_{t-1} %*% B_t  %*% t(1))
  
  ####
  # Set parameters
  n <- length(SL) # Length of time serie, refered as T in text
  t_1 <- matrix(1,nrow=2) # this is 1' or t(1)
  
  # Set parameters of state-dependent probabilities of observations
  # Uniform angle distribution for F: mu =0, kappa=0
  # Forward persistance for T: mu=0
  
  #####
  # Parameters to estimate
  
  # Note that to be able to use unconstrained optimizers
  # I have transformed the parameters
  # plogis constraint values betwen 0 & 1
  # exp() > 0 (actually, because of rounding >=0)
  gII <- plogis(x[1])
  gEE <- plogis(x[2])
  lI <- .Machine$double.xmin + exp(x[3])
  lE <- .Machine$double.xmin + exp(x[4])
  kE <- exp(x[5]) # Note that kappa can be 0 (uniform)
  
  
  ##
  # To allow to fix some of the parameters
  # This is need for the profile likelihood CI
  if(length(parF)>0){
    for (i in 1:length(parF)){
      assign(names(parF[i]),parF[[i]])
    }
  }
  
  # Gamma = the t.p.m
  # Only estimating gII, gEE,
  # a12 & a21 are (1-gII) and (1-gEE), respectively
  # This assures that the row of the t.p.m sum to 1
  # Which mean that you will end up in one of the 2 states at the next time step
  # with a probability of 1
  Gamma = matrix(c(gII,(1-gII),(1-gEE),gEE), nrow=2, byrow=T)
  
  # Here I am using the stationary distribution of the Markov Chain
  # for the initial distribution of the Markov Chain
  tpm <- matrix(c(gII, (1-gEE), (1-gII), gEE), nrow=2, byrow=TRUE)
  dI <- eigen(tpm)$vectors[,1]
  # Need to normalise the eigen vector to sum to 1
  dI <- dI/sum(dI)
  dI <- dI[1]
  delta <- c(dI, 1-dI)
  
  # P(x_t) is the diag of p_i(x_t)
  P <- function(SL_t,TA_t){
    # State 1: F
    # dexp(SL|lI) * dvonmises(TA|muI, kappa_F)
    pI <-   dexp(SL_t-SLmin,lI) * dvm(TA_t, 0, 0)
    # State 2: T
    # dexp(SL|lE) * dvonmises(TA|muE, kE)
    pE <-   dexp(SL_t-SLmin,lE) * dvm(TA_t, 0, kE)
    
    result <- diag(c(pI,pE))
  }
  
  
  ########################
  # Creating the variables that will save values during the recursive
  # calculation of the likelihood
  phi <- matrix(NA, ncol=2, nrow = n)
  
  #####################
  # Initialising the likelihood
  
  # The weight and scaled forward probabilities of the initial state
  w_1 <- delta %*% P(SL[1],TA[1]) %*% t_1
  phi[1,] <- (1/w_1) %*% delta %*% P(SL[1],TA[1])
  # the log likelihood of the 1st observation
  LL <- log(w_1)
  for (i in 2:n){
    v <- phi[(i-1),] %*% (Gamma %^% missL[i-1]) %*% P(SL[i],TA[i])
    u <- v %*% t_1
    LL <- LL + log(u)
    phi[i,] <- (1 / u) %*% v 
  }
  return(-LL) # Return loglikelihood
}


#######################################
# CCRW for numerical maximization with weibull and wrapped Cauchy

nllCCRWww <- function(SL,TA,x,parF=list(missL=missL)){
  # ParF should have missL
  
  ####
  # Set parameters
  n <- length(SL) # Length of time serie, refered as T in text
  t_1 <- matrix(1,nrow=2) # this is 1' or t(1)
  
  # Set parameters of state-dependent probabilities of observations
  # Uniform angle distribution for F: mu =0, kappa=0
  # Forward persistance for T: mu=0
  
  #####
  # Parameters to estimate
  
  # Note that to be able to use unconstrained optimizers
  # I have transformed the parameters
  # plogis constraint values betwen 0 & 1
  # exp() > 0 (actually, because of rounding >=0)
  gII <- plogis(x[1])
  gEE <- plogis(x[2])
  scI <- .Machine$double.xmin + exp(x[3])
  scE <- .Machine$double.xmin + exp(x[4])
  shI <- .Machine$double.xmin + exp(x[5])
  shE <- .Machine$double.xmin + exp(x[6])
  rE <- plogis(x[7]) # Note that kappa can be 0 (uniform)
  
  
  ##
  # To allow to fix some of the parameters
  # This is need for the profile likelihood CI
  if(length(parF)>0){
    for (i in 1:length(parF)){
      assign(names(parF[i]),parF[[i]])
    }
  }
  
  # Gamma = the t.p.m
  # Only estimating gII, gEE,
  # a12 & a21 are (1-gII) and (1-gEE), respectively
  # This assures that the row of the t.p.m sum to 1
  # Which mean that you will end up in one of the 2 states at the next time step
  # with a probability of 1
  Gamma = matrix(c(gII,(1-gII),(1-gEE),gEE), nrow=2, byrow=T)
  
  # Here I am using the stationary distribution of the Markov Chain
  # for the initial distribution of the Markov Chain
  tpm <- matrix(c(gII, (1-gEE), (1-gII), gEE), nrow=2, byrow=TRUE)
  dI <- eigen(tpm)$vectors[,1]
  # Need to normalise the eigen vector to sum to 1
  dI <- dI/sum(dI)
  dI <- dI[1]
  delta <- c(dI, 1-dI)
  
  # P(x_t) is the diag of p_i(x_t)
  P <- function(SL_t,TA_t){
    # Using weibull instead of exponential, shape (beta) and scale (eta) instead of lambda
    # Using wrapped Cauchy instead of vonMises,rho ineats of kappa. rho needs to be constrained between 0 & 1
    # The weibull distribution can only take values of 0 if the shI is 1, so remove SLmin, but still ignore the
    # SL that are 0
    
    # State 1: Intensive
    pI <-   dweibull(SL_t, shI, scI) * dwrpcauchy(TA_t, 0, 0)
    # State 2: Extensive
    pE <-   dweibull(SL_t, shE, scE) * dwrpcauchy(TA_t, 0, rE)
    
    result <- diag(c(pI,pE))
    return(result)
  }
  
  
  ########################
  # Creating the variables that will save values during the recursive
  # calculation of the likelihood
  phi <- matrix(NA, ncol=2, nrow = n)
  
  #####################
  # Initialising the likelihood
  
  # The weight and scaled forward probabilities of the initial state
  w_1 <- delta %*% P(SL[1],TA[1]) %*% t_1
  phi[1,] <- (1/w_1) %*% delta %*% P(SL[1],TA[1])
  # the log likelihood of the 1st observation
  LL <- log(w_1)
  for (i in 2:n){
    v <- phi[(i-1),] %*% (Gamma %^% missL[i-1]) %*% P(SL[i],TA[i])
    u <- v %*% t_1
    LL <- LL + log(u)
    phi[i,] <- (1 / u) %*% v 
  }
  return(-LL) # Return loglikelihood
}


#######################################
# LW
nllLW <- function(SL,TA,mu,parF=list('SLmin'=min(SL))){
	##
	# Parameters to estimate:
	# mu
	# Note that SLmin (a) is also considered an estimated parameter
	# According to Edwards (2008), to consider the full pdf of the measured movement
	# as potentially  power law, set a as the smallest measured movement

	##
	# To allow to fix some of the parameters
	# This is need for the profile/slice likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}
		
	##
	# Set parameters
	# Uniform turning angle distribution: mu=0, kappa=0

	LL <- log(mu-1) + (mu-1)*log(SLmin) - mu*log(SL) + log(dvm(TA, 0, 0))
	return(-sum(LL))
}

#######################################
# TLW

nllTLW <- function(SL,TA,mu,parF=list(SLmin=min(SL),SLmax=max(SL))){
	###
	# Parameter to estimate
	# mu
	# Note that SLmin (a) and SLmax (b) are also considered estimated parameters

	##
	# To allow to fix some of the parameters
	# This is need for the profile/slice likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}

	##
	# Set parameters
	# Uniform turning angle distribution: mu=0, kappa=0

	# from Edwards et al 2007 (supp): likelihood = (mu - 1) / (a^(1-mu) - b^(1-mu))^(-1) * SL^(-mu)
  # Unlike for the pure LW, the TLW can have mu <= 1, 
  # The likelihood for mu =1 is different than that for mu != 1, see Edwards et al. 2012
  # if mu = 1 
  if(mu == 1){
    LL <- log(1/(log(SLmax) - log(SLmin))) - mu*log(SL) + log(dvm(TA, 0, 0)) 
  }else{
    LL <- log((mu-1)/(SLmin^(1-mu) - SLmax^(1-mu))) - mu*log(SL) + log(dvm(TA, 0, 0)) 
  }
	return(-sum(LL))
}

#######################################
# BW
nllBW <- function(SL,TA,lambda,parF=list('SLmin'=min(SL))){
	##
	# Parameters to estimate:
	# lambda
	# Note that SLmin (a) is also considered an estimated parameter
	
	##
	# To allow to fix some of the parameters
	# This is need for the profile/slice likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}

	# Fixed parameters
	# mu <- 0 # mean of 0
	# kappa <- 0 # means that the distribution is uniform (no correlation in turning angle)

	LL <- dexp((SL-SLmin), lambda, log=TRUE) + log(dvm(TA, 0, 0))
	return(-sum(LL))
}

#######################################
# CRW
nllCRW <- function(SL,TA,kapp,lambda=parF$lambda,parF=list('SLmin'=min(SL))){
	##
	# Parameters to estimate:
	# lambda, kapp
	# Note that SLmin (a) is also considered an estimated parameter

	##
	# To allow to fix some of the parameters
	# This is need for the profile/slice likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}

	# Fixed parameters
	# For directional persistence: mu =0

	LL <- dexp((SL-SLmin), lambda, log=TRUE) + log(dvm(TA, 0, kapp))
	return(-sum(LL))
}

#######################################
# TBW

nllTBW <- function(SL,TA,x,parF=list('SLmin'=SLmin,'SLmax'=SLmax)){
	##
	# Parameters to estimate:
	# lambda
	# Note that SLmin (a) and SLmax (b) are also considered estimated parameters
	lambda <- exp(x)		
	# Fixed parameters
	# mu_B <- 0 # mean of 0
	# kappa_B <- 0 # means that the distribution is uniform (no correlation in turning angle)

	##
	# To allow to fix some of the parameters
	# This is need for the profile/slice likelihood CI
	if(length(parF)>0){
		for (i in 1:length(parF)){
			assign(names(parF[i]),parF[[i]])
		}
	}

	# According to Edwards 2011 supp
	# LL = n * log(lambda) - n*log(e^{-lambda*a} - e^{-lambda*b}) - lambda *sum(SL)
	LL <- (log(lambda) - log(exp(-lambda*SLmin)-exp(-lambda*SLmax)) - lambda*SL) + log(dvm(TA, 0,0))
	return(-sum(LL))
}

#######################################
# TCRW
nllTCRW <- function(SL,TA,lambda,kapp,SLmin=min(SL),SLmax=max(SL)){
	##
	# Parameters to estimate:
	# lambda, kapp
	# Note that SLmin (a) and SLmax are also considered an estimated parameter

	# Fixed parameter:
	# For directional persistence: mu= 0

	LL <- (log(lambda) - log(exp(-lambda*SLmin)-exp(-lambda*SLmax)) - lambda*SL) +
		log(dvm(TA, 0, kapp))
	return(-sum(LL))
}

######################
# CCRW HSMM (based on Langrock et al. 2012)

###
# function that derives the t.p.m. of the HMM that represents the HSMM
# dbinom reparametrise with mu
gen.Gamma.repar <- function(m,pSize,pMu){
  Gamma <- diag(m[1]+m[2])*0
  # p(r), for r=1,2,...N* (note that r=1 is 0 in dnbinom)
  #probs1 <- dnbinom(0:(m[1]-1),size=pSize[1],mu=pMu[1])
  #probs2 <- dnbinom(0:(m[2]-1),size=pSize[2],mu=pMu[2])
  probs1 <- dnbinom(0:(m[1]-1),size=pSize[1],prob=pMu[1])
  probs2 <- dnbinom(0:(m[2]-1),size=pSize[2],prob=pMu[2])
  
  # Denominator of c(r): 1 - sum_{k=1}^{r-1}p(k), so for r=1 -> 0 (because empty sum equals 0)
  # Use cumulative distribution function because it is the sum of the prob
  #den1 <- 1 - c(0,pnbinom(0:(m[1]-2),size=pSize[1],mu=pMu[1]))
  #den2 <- 1 - c(0,pnbinom(0:(m[2]-2),size=pSize[2],mu=pMu[2]))
  den1 <- 1 - c(0,pnbinom(0:(m[1]-2),size=pSize[1],prob=pMu[1]))
  den2 <- 1 - c(0,pnbinom(0:(m[2]-2),size=pSize[2],prob=pMu[2]))
  
  # To remove the chance of getting Inf
  probs1[which(den1<1e-12)] <- 1
  den1[which(den1<1e-12)] <- 1
  probs2[which(den2<1e-12)] <- 1
  den2[which(den2<1e-12)] <- 1
  
  # state aggregate 1
  Gamma[1:m[1],m[1]+1] <- probs1/den1 # c_1(r) for r=1,2,...,N_1* in first column of Beh 2
  diag(Gamma[1:(m[1]-1),2:m[1]]) <- 1-Gamma[1:(m[1]-1),m[1]+1] # 1-c_1(r), for r=1,2,...,N_1*-1
  Gamma[m[1],m[1]] <- 1 - Gamma[m[1],m[1]+1] # 1-c_1(N_1*)
  
  # state aggregate 2
  Gamma[m[1]+(1:m[2]),1] <- probs2/den2 # c_2(r) for r=1,2,...,N_2* in first column of Beh 1
  diag(Gamma[m[1]+1:(m[2]-1),m[1]+2:m[2]]) <- 1 - Gamma[m[1]+1:(m[2]-1),1] # 1-c_2(r), for r=1,2,...,N_2*-1
  Gamma[m[1]+m[2],m[1]+m[2]] <- 1 - Gamma[m[1]+m[2],1] # 1-c_2(N_2*)
  return(Gamma)
}


#########################################################
# Neg. log likelihood functions

nllHSMM <- function(SL, TA, x, parF){
  # oarf need missL and notMissLoc and m
  #####
  # Parameters to estimate - transforming for unconstrained parameter
  gSI <- .Machine$double.xmin + exp(x[1])
  gSE <- .Machine$double.xmin + exp(x[2])
#   gMI <- .Machine$double.xmin + exp(x[3])
#   gME <- .Machine$double.xmin + exp(x[4])
  gMI <- plogis(x[3])
  gME <- plogis(x[4])
  scI <- .Machine$double.xmin + exp(x[5])
  scE <- .Machine$double.xmin + exp(x[6])
  shI <- .Machine$double.xmin + exp(x[7])
  shE <- .Machine$double.xmin + exp(x[8])
  rE <- plogis(x[9]) # Note that kappa can be 0 (uniform)

  
  ##
  # To allow to fix some of the parameters
  # This is need for the profile likelihood CI
  if(length(parF)>0){
    for (i in 1:length(parF)){
      assign(names(parF[i]),parF[[i]])
    }
  }
  
  gamma <- gen.Gamma.repar(m,c(gSI,gSE),c(gMI,gME)) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  
  # This is the real length of the time series (include missing data)
  n <- length(SL) + sum(missL-1)

  obsProb <- matrix(rep(1,sum(m)*n),nrow=n)
  # Making a matrix with all the probability of observations for each state
  # Note that for missing location observation probablility is 1
  obsProb[notMisLoc,1:m[1]] <-  dweibull(SL, shI, scI) * dwrpcauchy(TA, 0, 0)
  obsProb[notMisLoc,(m[1]+1):sum(m)] <-  dweibull(SL, shE, scE) * dwrpcauchy(TA, 0, rE)

  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*obsProb[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    #if(any(is.nan(foo))){stop(print(i))}
  }
  return(-lscale)
}
