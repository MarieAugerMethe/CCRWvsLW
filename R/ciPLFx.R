#######################################
# Created By: Marie Auger-Methe
# Date created: November 21, 2011
# Updated: December, 2014

# This script contains the functions for computing the profile likelihood (PL) CI.
# For the PL CI you need the overall function that calculate the PL CI (CI.PL())
# for a given parameter,
# and the individual fx for each model that uses CI.PL()
# to calculate the PL CIs for all the parameters estimated through mle

#######################################
CI.PL.EM <- function(SL, TA, parI, parMLE, missL, parF, rang.b, mnll.m, notMisLoc, B=100, graph=TRUE){
  
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
    mnllRes <- EM_CCRW_HMM_CI(SL,TA,parMLE,missL,notMisLoc,parF=c(parF,rang.parI[i]))
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
  if(any(prof.lower > thr_95, na.rm=TRUE)){
    im <- max(which(prof.lower > thr_95))
    if(im >1){
      prof.lower[1:(im-1)] <- NA 
    }
  }
  if(any(prof.upper > thr_95, na.rm=TRUE)){
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
  
  if(graph==TRUE){
    # The likelihood profile is pretty ugly for the values at the boundary
    plot(rang.mnll~rang.parI, 
         main="", pch=20, ty="o",
         xlab=names(parI), ylab="Negative Log Likelihood")
    abline(h=thr_95, lty=2, col="darkgrey", lwd=2)
    abline(v=res, lty=3, col="darkgrey", lwd=2)
  }
  return(res)
}

#######################################
CI.PL <- function(SL, TA, parI, parMLE, trans.par, NLL, parF, rang.b, mnll.m, B=100, graph=TRUE, extOpt =FALSE){
  
  # Function inputs
  # SL: step length
  # TA: relative turn angle
  # parI: parameter of interest (for which the CI is computed with this fx) (not transformed)
  # parMLE.t: MLE values for all the parameters that are estimated through MLE (already transformed)
  # NLL: neg LL function for the model (one of the fx in nll.fx.R)
  # B: number of values of parI for which the min neg LL is estimated
  # rang.b: range boundaries for parI
  # mnll.m: min neg LL of the model when all parameters are minimised
  
  mnll_i <- floor(B/2)
  rang.parI <- numeric(B)
  rang.parI[1:mnll_i] <- seq(rang.b[1],parI,length=mnll_i)
  rang.parI[mnll_i:B] <- seq(parI,rang.b[2],length=(B-mnll_i+1))
  names(rang.parI) <- rep(names(parI),B)
  
  # Set memory for loop results
  rang.mnll <- numeric(B)
  
  for (i in 1:B){
    mnllRes <- tryCatch(optim(trans.par(parMLE),NLL,SL=SL,TA=TA, parF=c(parF,rang.parI[i])),
                        error=function(e) list('value'=NA, 
                                               'message'='Optim returned an error in the CI.PL function'))
    if(extOpt){
      if(!is.na(mnllRes$value)){
        mnllRes2 <- tryCatch(optim(mnllRes$par,NLL,SL=SL,TA=TA, parF=c(parF,rang.parI[i])),
                             error=function(e) list('value'=NA, 
                                                    'message'='Optim returned an error in the CI.PL function'))
        if(mnllRes2$value < mnllRes$value){mnllRes <- mnllRes2}
        mnllRes2 <- tryCatch(optim(mnllRes$par,NLL,SL=SL,TA=TA, parF=c(parF,rang.parI[i])),
                             error=function(e) list('value'=NA, 
                                                    'message'='Optim returned an error in the CI.PL function'))
        if(mnllRes2$value < mnllRes$value){mnllRes <- mnllRes2}
      }else{
        mnllRes <- tryCatch(optim(trans.par(parMLE)*1.1,NLL,SL=SL,TA=TA, parF=c(parF,rang.parI[i])),
                            error=function(e) list('value'=NA, 
                                                   'message'='Optim returned an error in the CI.PL function'))
      }
    }
    
    rang.mnll[i] <- mnllRes$value
    mnllRes$message
  }
  
  # Step 2. Divide into two branch (Lower and upper)
  prof.lower <- rang.mnll[1:mnll_i]
  prof.l.avec <- rang.parI[1:mnll_i]
  prof.upper <- rang.mnll[mnll_i:B]
  prof.u.avec <- rang.parI[mnll_i:B]
  
  # Step 3. Using linear interpolation to get the value of X for which
  # neg LL (for the parameter range) = neg LL (from the real MLE) + chisqdist(0.95)/2
  thr_95 <- mnll.m + qchisq(0.95,1)/2
  
  # So we get to good peak
  if(any(prof.lower > thr_95, na.rm=TRUE)){
    im <- max(which(prof.lower > thr_95))
    if(im >1){
      prof.lower[1:(im-1)] <- NA 
    }
  }
  if(any(prof.upper > thr_95, na.rm=TRUE)){
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
  
  if(graph==TRUE){
    # The likelihood profile is pretty ugly for the values at the boundary
    plot(rang.mnll~rang.parI, 
         main="", pch=20, ty="o",
         xlab=names(parI), ylab="Negative Log Likelihood")
    abline(h=thr_95, lty=2, col="darkgrey", lwd=2)
    abline(v=res, lty=3, col="darkgrey", lwd=2)
  }
  return(res)
}

#######################################
# Becaus ethe other models are not using an EM algorthim and can generally be divided into
# indepedent parameters you can do a slice.
CI.Slice <- function(SL, TA, parI, parF, trans.par, NLL, rang.b, mnll.m, B=100, graph=TRUE){
  
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
  
  # Creating a matrix that will save the different minimization results
  rang.mnll <- numeric(B)
  
  for (i in 1:B){
    rang.mnll[i] <- NLL(SL=SL,TA=TA,trans.par(rang.parI[i]),parF=parF)
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
  
  if(graph==TRUE){
    # The likelihood profile is pretty ugly for the values at the boundary
    plot(rang.mnll~rang.parI,
         main="",pch=20, ty="o",
         xlab=names(parI), ylab="Negative Log Likelihood")
    abline(h=(mnll.m + qchisq(0.95,1)/2), lty=2, col="darkgrey", lwd=2)
    abline(v=res, lty=3, col="darkgrey", lwd=2)
  }
  
  return(res)
}

#######################################
# CCRW_HMM

ciCCRWpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  parF <- list('SLmin'= movD$SLmin)
  
  CI <- matrix(NA,nrow=5,ncol=3)
  rownames(CI) <- names(mleM[1:5])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  CI[,1] <- mleM[1:5]
  
  if(graph==TRUE){
    layout(matrix(1:5, nrow=1))
  }
  
  for(i in 1:5){
    CI[i,2:3] <- CI.PL.EM(movD$SL, movD$TA, mleM[i], mleM[1:6], movD$missL, parF, rangePar[i,], mleM['mnll'], 
                          movD$notMisLoc, B, graph)  
  }
  
  return(CI)
}

#######################################
# CCRW - numerical maximization

ciCCRWdnpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  parF <- list('SLmin'= movD$SLmin, 'missL' = movD$missL)

  CI <- matrix(NA,nrow=5,ncol=3)
  rownames(CI) <- names(mleM[1:5])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  CI[,1] <- mleM[1:5]
  
  if(graph==TRUE){
    layout(matrix(1:5, nrow=1))
  }
  
  for(i in 1:5){
    CI[i,2:3] <- CI.PL(movD$SL, movD$TA, mleM[i], mleM[1:5], transParCCRWdn, nllCCRWdn, parF, rangePar[i,],
                       mleM['mnll'], B=B, graph)
  }
  return(CI)
}

ciCCRWwwpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  parF <- list('missL' = movD$missL)
  
  CI <- matrix(NA,nrow=7,ncol=3)
  rownames(CI) <- names(mleM[1:7])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  CI[,1] <- mleM[1:7]
  
  if(graph==TRUE){
    layout(matrix(1:7, nrow=1))
  }
  
  for(i in 1:7){
    CI[i,2:3] <- CI.PL(movD$SL, movD$TA, mleM[i], mleM[1:7], transParCCRWww, nllCCRWww, parF, rangePar[i,],
                       mleM['mnll'], B=B, graph)  
  }
  
  return(CI)
}

ciHSMMpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  parF <- list("missL"= movD$missL, "notMisLoc"= movD$notMisLoc, "m"=c(10,10))
  
  CI <- matrix(NA,nrow=9,ncol=3)
  rownames(CI) <- names(mleM[1:9])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  CI[,1] <- mleM[1:9]
  
  if(graph==TRUE){
    layout(matrix(1:9, nrow=1))
  }
  
  for(i in 1:9){
    CI[i,2:3] <- CI.PL(movD$SL, movD$TA, mleM[i], mleM[1:9], transParHSMM, nllHSMM, parF, rangePar[i,],
                       mleM['mnll'], B=B, graph, extOpt = TRUE) 
  }
  
  return(CI)
}

# HSMM with gPI=gPE

ciHSMMlpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0, nPar=8){
  movD <- movFormat(movltraj, TAc)
  parF <- list("missL"= movD$missL, "notMisLoc"= movD$notMisLoc, "m"=c(10,10))
  
  CI <- matrix(NA,nrow=nPar,ncol=3)
  rownames(CI) <- names(mleM[1:nPar])
  colnames(CI) <- c("estimate","L95CI", "U95CI")
  CI[,1] <- mleM[1:nPar]
  
  if(graph==TRUE){
    layout(matrix(1:nPar, nrow=1))
  }
  
  for(i in 1:nPar){
    CI[i,2:3] <- CI.PL(movD$SL, movD$TA, mleM[i], mleM[1:nPar], transParHSMMl, nllHSMMl, parF, rangePar[i,],
                       mleM['mnll'], B=B, graph, extOpt = TRUE) 
  }
  
  return(CI)
}



#######################################
# LW

ciLWpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL  
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  
  parF <- list('SLmin'=SLmin)
  # Some models used transformed parameters in nll
  # Thus trans.par called in an input of CI.slice
  trans.par <- function(x){x}
  
  # The PL likelihood is more a slice in this case
  # since the other parameter estimated: a
  # is fixed to its value
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI_PL", "U95CI_PL")
  
  CI[,1] <- mleM[1]
  CI[,2:3] <- CI.Slice(SL, TA, mleM[1], parF, trans.par, nllLW, rangePar, mleM['mnll'], B, graph)
  
  return(CI)
}

#######################################
# TLW

ciTLWpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL  
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  
  parF <- list('SLmin'=SLmin,'SLmax'=SLmax)
  
  # The PL likelihood is more a slice in this case
  # since the other parameter estimated: a
  # is fixed to its value
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI_PL", "U95CI_PL")
  
  CI[,1] <- mleM[1]
  CI[,2:3] <- CI.Slice(SL, TA, mleM[1], parF, transParx, nllTLW, rangePar, mleM['mnll'], B, graph)
  
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

ciEpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL  
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  
  parF <- list('SLmin'=SLmin)
  
  # The PL likelihood is more a slice in this case
  # since the other parameter estimated: a
  # is fixed to its value
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI_PL", "U95CI_PL")
  
  CI[,1] <- mleM[1]
  CI[,2:3] <- CI.Slice(SL, TA, mleM[1], parF, transParx, nllBW, rangePar, mleM['mnll'], B, graph)
  
  return(CI)
}

#######################################
# K
# This is actually based on nll.CRW

ciKpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL  
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  
  parF <- list('SLmin'=SLmin, 'lambda'=mleM[1])

  # The PL likelihood is more a slice in this case
  # since the other parameter estimated: a
  # is fixed to its value
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[2])
  colnames(CI) <- c("estimate","L95CI_PL", "U95CI_PL")
  
  CI[,1] <- mleM[2]
  CI[,2:3] <- CI.Slice(SL, TA, mleM[2], parF, transParx, nllCRW, rangePar, mleM['mnll'], B, graph)
  
  return(CI)
}

#######################################
# TE
# This is actually based on nll.TBW

ciTEpl <- function(movltraj, mleM, rangePar, B=100, graph=TRUE, TAc=0){
  movD <- movFormat(movltraj, TAc)
  SL <- movD$SL
  TA_C <- movD$TA_C
  TA <- movD$TA
  SLmin <- movD$SLmin
  SLmax <- movD$SLmax
  missL <- movD$missL  
  n <- movD$n
  notMisLoc <- movD$notMisLoc
  
  parF <- list('SLmin'=SLmin,'SLmax'=SLmax)
  
  CI <- matrix(NA, nrow=1, ncol=3)
  rownames(CI) <- names(mleM[1])
  colnames(CI) <- c("estimate","L95CI_PL", "U95CI_PL")
  
  CI[,1] <- mleM[1]
  CI[,2:3] <- CI.Slice(SL, TA, mleM[1], parF, transParTBW, nllTBW, rangePar, mleM['mnll'], B, graph)  

  return(CI)
}
