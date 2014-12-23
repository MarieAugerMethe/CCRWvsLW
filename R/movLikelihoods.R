#######################################
# Created By: Marie Auger-Methe
# Date created: March 28, 2011
# Updated: April 3, 2013

movLikelihoods <- function(movltraj, graph=TRUE, PRdetails=FALSE, TAc=0){
  
  #######################################
  # This script estimates the parameters and calculates the AIC of multiple movement models:
  # 1. CCRW: combined BW and CRW for which the states are following a Markov process 
  # 2. LW: Levy walk
  # 3. TLW: Truncated Levy walk
  # 4. BW: Brownian walk
  # 5. CRW: Correlated random walk
  # 6. TBW: Truncated Brownian walk
  # 7. TCRW: Truncated Correlate random walk
  
  # It uses maximum likelihood paradigm for parameter estimation
  # For the likelihood function we are using both the angle and the step length distributions.
  # Because it helps differentiate between the 2 states and because it has been recommended 
  # to incorporate the turning in the analysis (Plank et al. 2011)
  
  # The angle is the relative turning angle. 
  # This allows me to model correlation in turning angle 
  # without violating the independence assumption of the likelihood
  # The relative turning angle models directionality by being
  # peaked at zero.
  
  # The CCRW assumes that the "foraging behavior" is a Brownian (non-correlated) walk
  # and that "travelling behaviour" has some correlation in turning angle
  
  # Because the likelihood function of the CCRW now includes the turning angle
  # and AIC comparison requires the model to use exactly the same data
  # I have also re-wrote the likelihood function for the LW (and TLW)
  # to include the turning angle distribution.
  # The likelihood functions for LW and TLW orignally only applied to the step length distribution.
  # The LW and TLW are modeled with a uniform distribution for the turning angle
  
  # To be consitent with previous studies
  # and because we are eliminating the steps at the exact same location
  # we use a minimum step length for all models.
  
  # Following Ben Bolker's book and Zucchini and MacDonald 2009:
  # 1. I have used the nlm optimization method for functions with multiple parameters to estimate (what Z&M 2009 use for HMM)
  #  I have chose nlm() instead of optim because it is faster
  # 	nlm() will return warning messages everytime an NA/Inf is encountered, but if there is not too many of these
  #	it shouldn't affect the minimization
  # 2. I have use optimize for the function with only one parameter to estimate
  # 3. To limit the parameter to values that are sensical 
  #    I have transformed the parameters (either by exponentiated or logit transform )
  
  #######################################
  # All the packages and functions required are loaded with sourcingFx.R
  
  #######################################
  # Formating movement trajectory
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL	
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  movltraj <- movD$movltraj
  
  ######################
  # Finding the mle
  
  mleCCRW <- mnllCCRW(SL, TA, TA_C, missL, notMisLoc, SLmin)
  
  mleLW <- mnllLW(SL, TA, SLmin)
  
  mleTLW <- mnllTLW(SL, TA, SLmin, SLmax)
  
  mleBW <- mnllBW(SL, TA, SLmin)
  
  mleTBW <- mnllTBW(SL, TA, SLmin, SLmax, mleBW['lambda'])
  
  mleCRW <- mnllCRW(SL, TA_C, TA, SLmin, mleBW['lambda'])
  
  # Lambda is the same as lambda_TBW
  # and kappa is the same as kappa_CRW
  mleTCRW <- mnllTCRW(SL, TA_C, TA, SLmin, SLmax, mleTBW['lambda'], mleCRW['kappa'])
  
  # Table with overall results from mle
  mleMov <- list(mleCCRW, mleLW, mleTLW, mleBW, mleCRW, mleTBW, mleTCRW)
  names(mleMov) <- c("CCRW", "LW", "TLW", "BW", "CRW", "TBW", "TCRW")

  ######################
  # Calculating the CI
  
  # CI for each of the parameter that were estimated using the mle.
  # Note that I'm not getting a CI for 'a' and 'b'
  # I'm calculating quadractic approximation of the CI 
  # using the Hessian or the second derivative of the likelihood at the mle.
  
  # Creating a list for the CI
  CI <- vector('list', length=6)
  names(CI) <- c(names(mleMov[1:3]),"E", "TE", "K")
  
  # CCRW
  CI[['CCRW']] <- ciCCRW(SL, TA, SLmin, missL, notMisLoc, mleCCRW)
  
  # LW
  CI[['LW']] <- ciLW(SL, TA, SLmin, mleLW)
  
  # TLW
  CI[['TLW']] <- ciTLW(SL,TA, SLmin, SLmax, mleTLW)
  
  # E
  CI[['E']] <- ciE(SL, TA, SLmin, mleBW)
  
  # TE
  CI[['TE']] <- ciTE(SL, TA, SLmin, SLmax, mleTBW)
  
  # K
  CI[['K']] <- ciK(SL, TA, SLmin, mleCRW)
  
  #######
  # Test of absolute fit
	pseudoRes <- pseudo(SL, TA_C, TA, SLmin, SLmax, missL, notMisLoc, n,
                      mleMov, PRdetails, graph)
  
	if(graph==TRUE){
    #windows()
    layout(matrix(1:3,nrow=1))
    # Movemeth trajectory with CCRW 
    gamm <- matrix(c(mleMov$CCRW[1], 1-mleMov$CCRW[2], 1-mleMov$CCRW[1], mleMov$CCRW[2]),2)
    w <- HMMwi(SL,TA,missL,SLmin, mleMov$CCRW[3:4], mleMov$CCRW[5], gamm,
                mleMov$CCRW[6:7], notMisLoc)
    plot(movltraj, addpoints=FALSE, final=FALSE)
    points(movltraj[[1]]$x,movltraj[[1]]$y,bg=c(0,grey(w[2,]),0),pch=21, cex=c(0,w[1,],0)+0.7)
		# histogram with fit
    fitGraph(SL, TA, SLmin, SLmax, n, mleMov)
	}

	res <- list(mleMov, pseudoRes, CI)
	names(res) <- c("mleMov","pseudoRes", "CI")
	return(res)

}
