#######################################
# Created By: Marie Auger-Methe
# Date created: March 28, 2011
# Updated: April 3, 2013

movLikelihoods <- function(movltraj, graph=TRUE, PRdetails=FALSE, TAc=0, conts=TRUE, dn=FALSE, 
                           ww=FALSE, hs=FALSE, hsl=FALSE, hss=FALSE, hsp=FALSE, hspo=FALSE){
  
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
  if(dn){
    mleCCRW <- mnllCCRWdn(SL, TA, TA_C, missL, SLmin)
  }else{
    mleCCRW <- mnllCCRW(SL, TA, TA_C, missL, notMisLoc, SLmin) 
  }
  
  mleLW <- mnllLW(SL, TA, SLmin)
  
  mleTLW <- mnllTLW(SL, TA, SLmin, SLmax, conts)
  
  mleBW <- mnllBW(SL, TA, SLmin)
  
  mleTBW <- mnllTBW(SL, TA, SLmin, SLmax, mleBW['lambda'])
  
  mleCRW <- mnllCRW(SL, TA_C, TA, SLmin, mleBW['lambda'])
  
  # Lambda is the same as lambda_TBW
  # and kappa is the same as kappa_CRW
  mleTCRW <- mnllTCRW(SL, TA_C, TA, SLmin, SLmax, mleTBW['lambda'], mleCRW['kappa'])
  
  # Table with overall results from mle
  mleMov <- list(mleCCRW, mleLW, mleTLW, mleBW, mleCRW, mleTBW, mleTCRW)
  names(mleMov) <- c("CCRW", "LW", "TLW", "BW", "CRW", "TBW", "TCRW")

  if(ww){
    mleCCRWww <- mnllCCRWww(SL, TA, TA_C, missL)
    mleMov$CCRWww <- mleCCRWww
  }
  
  if(hs){
    mleHSMM <- mnllHSMM(SL, TA, TA_C, missL, notMisLoc)  
    mleMov$HSMM <- mleHSMM
  }
  
  if(hsl){
    mleHSMMl <- mnllHSMMl(SL, TA, TA_C, missL, notMisLoc)  
    mleMov$HSMMl <- mleHSMMl
  }
  
  if(hss){
    mleHSMMs <- mnllHSMMs(SL, TA, TA_C, missL, notMisLoc)  
    mleMov$HSMMs <- mleHSMMs
  }
  
  if(hsp){
    mleHSMMp <- mnllHSMMp(SL, TA, TA_C, missL, notMisLoc)  
    mleMov$HSMMp <- mleHSMMp
  }
  
  if(hspo){
    mleHSMMpo <- mnllHSMMpo(SL, TA, TA_C, missL, notMisLoc)  
    mleMov$HSMMpo <- mleHSMMpo
  }
  
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
  if(dn){
    CI[['CCRW']] <- ciCCRWdn(SL, TA, SLmin, missL, mleCCRW)
  }else{
    CI[['CCRW']] <- ciCCRW(SL, TA, SLmin, missL, notMisLoc, mleCCRW)
  }
  
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
  
  if(ww){
    CI$CCRWww <- ciCCRWww(SL,TA,missL,mleCCRWww)
  }
  
  if(hs){
    CI$HSMM <- ciHSMMg(SL, TA, missL, notMisLoc,
                       mleHSMM, 9, nllHSMM, transParHSMM)
  }
  
  if(hsl){
    CI$HSMMl <- ciHSMMg(SL, TA, missL, notMisLoc,
                        mleHSMMl, 8, nllHSMMl, transParHSMMl)
  }
  
  if(hss){
    CI$HSMMs <- ciHSMMg(SL, TA, missL, notMisLoc,
                        mleHSMMs, 8, nllHSMMs, transParHSMMs)
  }
  
  if(hsp){
    CI$HSMMp <- ciHSMMg(SL, TA, missL, notMisLoc,
                        mleHSMMp, 7, nllHSMMp, transParHSMMp)
  }
  
  if(hspo){
    CI$HSMMpo <- ciHSMMg(SL, TA, missL, notMisLoc,
                        mleHSMMpo, 5, nllHSMMpo, transParHSMMpo)
  }
  
  #######
  # Test of absolute fit
	pseudoRes <- pseudo(SL, TA_C, TA, SLmin, SLmax, missL, notMisLoc, n,
                      mleMov, PRdetails, graph, dn=dn, ww=ww, hs=hs, hsl=hsl, hsp=hsp, hspo=hspo)
  
  
  
	if(graph==TRUE){
	  aiccAllccrw <- matrix(NA, ncol=length(CI)-5)
	  colnames(aiccAllccrw) <- c("o", "dn", "ww", "hs", "hsl", "hss", "hsp", "hspo")[c(!dn,dn,ww,hs,hsl,hss,hsp,hspo)]
	  aiccAllccrw[1] <- mleMov[[1]]["AICc"]
	  if(ncol(aiccAllccrw)>1){
	    for(i in 2:ncol(aiccAllccrw)){
	      aiccAllccrw[i] <- mleMov[[i+6]]["AICc"]
	    }  
	  }
	  
	  # Best ccrw
	  ccrwBm <- colnames(aiccAllccrw)[which.min(aiccAllccrw)]
	  
# 	  # Plot the best CCRW if ww explored
# 	  if(ww){
# 	    wwB <- which.min(c(mleMov$CCRW["AICc"], mleMov$CCRWww["AICc"])) == 2
# 	  }else{
# 	    wwB <- FALSE
# 	  }
    
    if(ccrwBm == "o"){
      gamm <- matrix(c(mleMov$CCRW[1], 1-mleMov$CCRW[2], 1-mleMov$CCRW[1], mleMov$CCRW[2]),2)
      w <- HMMwi(SL,TA,missL,SLmin, mleMov$CCRW[3:4], mleMov$CCRW[5], gamm,
                 mleMov$CCRW[6:7], notMisLoc)
    }else if(ccrwBm == "dn"){
      gamm <- matrix(c(mleMov$CCRW[1], 1-mleMov$CCRW[2], 1-mleMov$CCRW[1], mleMov$CCRW[2]),2)
      w <- HMMwi(SL,TA,missL,SLmin, mleMov$CCRW[3:4], mleMov$CCRW[5], gamm,
                 mleMov$CCRW[8:9], notMisLoc)
    }else if(ccrwBm == "ww"){
      gamm <- matrix(c(mleMov$CCRWww[1], 1-mleMov$CCRWww[2], 1-mleMov$CCRWww[1], mleMov$CCRWww[2]),2)
      w <- HMMwiww(SL, TA, missL, 
                   mleMov$CCRWww[5:6], mleMov$CCRWww[3:4], mleMov$CCRWww[7], gamm, mleMov$CCRWww[9:10],
                   notMisLoc)
    }else if(ccrwBm == "hs"){
      w <- HSMMwi(SL, TA, missL, notMisLoc, mleMov$HSMM[1:2], mleMov$HSMM[3:4], mleMov$HSMM[5:6], mleMov$HSMM[7:8], mleMov$HSMM[9])
    }else if(ccrwBm == "hsl"){
      w <- HSMMwi(SL, TA, missL, notMisLoc, mleMov$HSMMl[1:2], mleMov$HSMMl[rep(3,2)], mleMov$HSMMl[4:5], mleMov$HSMMl[6:7], mleMov$HSMMl[8])
    }else if(ccrwBm == "hsp"){
      w <- tryCatch(HSMMpwi(SL, TA, missL, notMisLoc, mleMov$HSMMp[1:2], mleMov$HSMMp[3:4], mleMov$HSMMp[5:6], mleMov$HSMMp[7]),
                      error=function(e){matrix(0, nrow=2,ncol=notMisLoc[length(notMisLoc)])})
    }else if(ccrwBm == "hspo"){
      w <- tryCatch(HSMMpowi(SL, TA, missL, notMisLoc, mleMov$HSMMpo[1:2], mleMov$HSMMpo[3:4], mleMov$HSMMpo[5]),
                       error=function(e){matrix(0, nrow=2,ncol=notMisLoc[length(notMisLoc)])})
    }
    
	  layout(matrix(1:3,nrow=1))    
    plot(movltraj, addpoints=FALSE, final=FALSE)
    points(movltraj[[1]]$x,movltraj[[1]]$y,bg=c(0,grey(w[2,]),0),pch=21, cex=c(0,w[1,],0)+0.7)
		# histogram with fit
    fitGraph(SL, TA, SLmin, SLmax, n, mleMov, ccrwBm)
	}

	res <- list(mleMov, pseudoRes, CI)
	names(res) <- c("mleMov","pseudoRes", "CI")
	return(res)

}
