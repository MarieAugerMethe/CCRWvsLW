# I'm using the normal pseudo residuals to test the absolute fit of the model
# I'm getting the pseudo residuals of the SL and the TA independently
# and this is because I'm assuming that SL and TA are independent
# This assumption is not quite true for the CCRW
# But I still think it's ok because th dependence is mainly
# associaed with the behaviour and in both case the behaviour
# is estimated with both the SL and TA
# There are some special tricks with the CCRW see below for more info

# I'm using Zucchini and MAcDonald 2009 definition of pseudoresiduals
# The general idea is:
# 1. I'm getting the uniform pseudo residuals
# Which are simply the cummulative density function value for observation x_t
# u_t = F(x_t) = Pr[X_t<=x_t]
# in the case of CCRW it is more complicated
# they are defined as u_t = F(x_t) = Pr[X_t<=x_t | X^{-t}=x^{-t}]
# So it is based on the conditional distribution

# 2. I'm transforming the unfirom pseudo residuals into normal pseudo residual
# which is simply using inverse function of a standard normal (or quantile fx)
# z_t = phi{-1}(u_t)

# 3. I'm using a shapiro test of normality to test whether the normal 
# pseudo-residuals are normal

################################################################
# Functions needed

###
# test of uniformity based on G-test

pseudo.u.test <- function(U, nP, n, graph, movMes){

  # Binning the uniform pseudo residuals
  # Dividing the uniform distribution in as many bins as possible
  # that will have 5 expected values.
  # Since the degrees of freedom depends on the number
  # of bins and the number of parameters estimated
  # df = nB - nP -1
  # and the biggest number of parameter estimated is 3
  # we need a minium of 5 bins with 10 expected
  # so a minimum of 50 data points.
  
  # Constant size bin
  nB <- floor(n/10)
  U_E <- n/nB
  U_B <- seq(0,1,length.out=(nB+1))
  U_B <- hist(U,breaks=U_B,plot=FALSE)$counts
  
  ##
  # The actual G-test
  # According to Sokal and Rohlf (1995)
  # GVal <- 2*sum(O*log(O/E))
  GVal <- U_B*log(U_B/U_E)
  GVal[U_B==0] <- 0
  GVal <- 2*sum(GVal)
  
  # The df depends only on the non-empty bins
  nB <- nB - sum(U_B==0)
  
  # According to Sokal and Rohlf 1995
  # the degrees of freedom depends on whether 
  # there is an intrinsic hypothesis
  # (a hypothesis based on the data) or not
  # But since I'm doing the analysis on the pseudo residuals
  # I'm not sure whether I should penalise
  # the test for the number of parameters
  dfGT <- nB - nP - 1
  
  # Williams correction according to Sokal and Rohlf (1995)
  # q=1+(a^2-1)/6* n*dfGT, where a is number of category
  wc <- 1 + (nB^2 -1)/(6 * n * dfGT)
  
  pValGT <- pchisq(GVal/wc,dfGT,lower.tail=FALSE)
  pValGT <- rbind(Gcor=(GVal)/wc,pValGT)
  
  if(graph==TRUE){
    par(mfrow=c(1,2))
    hist(U, nB, xlab = "Uniform pseudo-residuals", ylab=paste("Frequency",movMes), main="")
    box()
    acf(U, xlab="Lag (steps)", ylab=paste("ACF",movMes), main="")
  }
  return(pValGT)
}

###
# Function to combine p-values from the 2 test of fit (SL & TA)

# According to Quinn & Keough, you can combined p-value
# with comb_stat = -2* sum_{i=1}^{c}(ln(p_c)),
# where p_c is the p value of the test c
# and comb_stat is chi squared distributed with 2*c degrees of freedom
# This is often refered as the Fisher's method

combP <- function(p1,p2){
	h <- -2*(log(p1)+log(p2))
	pc <- pchisq(h,df=4,lower.tail=F)
	res <- cbind(h,pc)
	colnames(res) <- c('combStat','pval')
	return(res)
}

###
# For most function there is a function to get the cummulative
# density function but for some of them I had to write the function
##

# Probability density function of TE
dtexp <- function(SL,SLmin,SLmax,parTE){
  parTE * exp(-parTE *SL)/ (exp(-parTE*SLmin) - exp(-parTE*SLmax))
}

ptexp <- function(SL,SLmin,SLmax,parTE){

  sapply(SL,
         function(x) integrate(dtexp,SLmin,x,SLmin=SLmin,SLmax=SLmax,parTE=parTE)[[1]])
}



##
# CCRW

# I'm using the ordinary pseudo residuals for CCRW defined by Zucchini and MacDonald 2009
# The only difference is that I'm using both SL and TA
# to get the w associated with the conditional probabilities
# but only p(x_t) appropriate to the good measure (SL or TA)
# when I'm calculating the pseudo-residuals


# function that test the absolut fit

pseudo <- function(SL,TA_C,TA,SLmin,SLmax,missL,notMisLoc,n,mleMov,PRdetails,graph,Uinfo=FALSE){
  
	##########################################################
	# Pseudo residuals

  # If graph=TRUE get the best model
  # since only the pseudo residuals of the best model are shown
  if(graph){
    bestM <- grep("AICc",names(unlist(mleMov)))
    # Order: CCRW, LW, TLW, BW, CRW, TBW, TCRW
    bestM <- which.min(unlist(mleMov)[bestM])
    graph <- rep(FALSE,7)
    graph[bestM] <- TRUE
  }else{
    graph <- rep(FALSE,7)
  }
  
  # By definition the SL that are corresponding to SLmin and SLmax
  # will be outliers because it is very unlikely o have values that are
  # exactly at the truncation value
  # Because the truncation values are set to be the smallest and greatest observed
  # It cause issues with normality test,
  # It is because the normal pseudo residuals of SLmin and SLmax
  # will be -Inf and Inf
  # SO the easiest way around this is to remove these to outliers
  
  Z <- matrix(NA,2,8)
  colnames(Z) <- c('SL_CCRW', 'SL_LW', 'SL_TLW','SL_E','SL_TE',
                   'TA_CCRW', 'TA_U','TA_VM')
  rownames(Z) <- c('W','pval')
  
	##
	# CCRW
  
  if(!any(is.na(mleMov$CCRW[1:6]))){
    gamm <- matrix(c(mleMov$CCRW[1], 1-mleMov$CCRW[2], 1-mleMov$CCRW[1], mleMov$CCRW[2]),2)
    w <- HMMwi(SL,TA,missL,SLmin,
      mleMov$CCRW[3:4], mleMov$CCRW[5], gamm, mleMov$CCRW[6:7], notMisLoc)
    
    U_SL_CCRW <- w[1,notMisLoc] * pexp(SL-SLmin,mleMov$CCRW[3]) +
      w[2,notMisLoc] * pexp(SL-SLmin,mleMov$CCRW[4])
    
    U_TA_CCRW <- w[1,notMisLoc] * pvonmises(TA_C, circular(0), 0) +
      w[2,notMisLoc] * pvonmises(TA_C, circular(0), mleMov$CCRW[5])
    
    # Although I'm estimated 7 parameters in total for the CCRW
    # and all parameters are used to identify the expected behavioral state
    # I'm only taking into acount the parameters that trully are associated with
    # each set of pseudo-residuals.
    # SO for SL, you have 3 parameters: minSL, lambda_T, lambda_F
    # SO for TA, you have 1 parameter: kappa_T
    Z[,1] <- pseudo.u.test(U_SL_CCRW, 3, n, graph[1], "SL")
    Z[,6] <- pseudo.u.test(U_TA_CCRW, 1, n, graph[1], "TA")
  }

  
  #####
  # SL
  
	##
	# LW

	U_LW <- ppareto(SL,SLmin,mleMov$LW[1]-1)

	##
	# TLW

	U_TLW <- ptruncpareto(SL,SLmin,SLmax,mleMov$TLW[1]-1) # Used to be ptpareto

	##
	# E
	# For both BW and CRW

	U_E <- pexp(SL-SLmin,mleMov$BW[1])

	##
	# TE
	# For both TBW and TCRW

	U_TE <- ptexp(SL,SLmin,SLmax,mleMov$TBW[1])


	#####
	# TA
  
	##
	# U
	# Circular uniform
	# For LW, TLW, BW, and TBW

  U_U <- pvonmises(TA_C,circular(0),0)
  
	##
	# VM
	# von Mises
	# For CRW and TCRW

  U_VM <- pvonmises(TA_C,circular(0),mleMov$CRW[2])
   
  
	##########################################################
	# Normal pseudo residuals



	######################
	# SL
	
  # nP = SLmin + mu
	Z[,2] <- pseudo.u.test(U_LW,2,n,graph[2],"SL")
  # nP = SLmin + mu + SLmax
	Z[,3] <- pseudo.u.test(U_TLW,3,n,graph[3],"SL")
  # nP = SLmin + lambda
	Z[,4] <- pseudo.u.test(U_E,2,n,graph[4]|graph[5],"SL")
  # nP = SLmin + lambda + SLmax
	Z[,5] <- pseudo.u.test(U_TE,3,n,graph[6]|graph[7],"SL")

	######################
	# TA
  
  # nP = nothing!
	Z[,7] <- pseudo.u.test(U_U,0,n,graph[2]|graph[3]|graph[4]|graph[6],"TA")
  # nP = kapa
	Z[,8] <- pseudo.u.test(U_VM,1,n,graph[5]|graph[7],"TA")

	PR <- matrix(NA,2,7)
	colnames(PR) <- c('CCRW', 'LW', 'TLW', 'BW', 'CRW', 'TBW', 'TCRW')
	rownames(PR) <- c('CombStat','pval')
	PR[,1] <- combP(Z[2,1],Z[2,6]) # CCRW
	PR[,2] <- combP(Z[2,2],Z[2,7]) # LW
	PR[,3] <- combP(Z[2,3],Z[2,7]) # TLW
	PR[,4] <- combP(Z[2,4],Z[2,7]) # BW
	PR[,5] <- combP(Z[2,4],Z[2,8]) # CRW
	PR[,6] <- combP(Z[2,5],Z[2,7]) # TBW
	PR[,7] <- combP(Z[2,5],Z[2,8]) # TCRW

  if(PRdetails & Uinfo){
    return(list('PR'=PR,'Z'=Z, 
                'U'= cbind(U_SL_CCRW, U_LW, U_TLW, U_E, U_TE,
                           U_TA_CCRW, U_U, U_VM)))
  }else if(PRdetails){
    return(list('PR'=PR,'Z'=Z))
  }else if(Uinfo){
    return(list('PR'=PR, 
                'U'= cbind(U_SL_CCRW, U_LW, U_TLW, U_E, U_TE,
                           U_TA_CCRW, U_U, U_VM)))
  }else{
    return(PR) 
  }
	
}