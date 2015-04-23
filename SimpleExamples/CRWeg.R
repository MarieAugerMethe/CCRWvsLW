library(CCRWvsLW)
l <- 0.01
k <- 10
a <- 1 
mov <- simmCRW(500,l,k,a)
movRes <- movLikelihoods(mov, PRdetails=TRUE)
movRes

# To look at the best model compare AICc
AICcRes <- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] # CRW is the best
AICcRes - min(AICcRes) # By far, but CCRW is also good (normal since it can reproduce the CRW)

# Compare resuts to simulation values
cbind(movRes$mle$CRW[1:3], c(l,k,a)) # Pretty good

# Look at test of absolute fit
# With an alpha of 0.05, not significantly different from CCRW, CRW, TCRW
round(movRes$pseudoRes$PR["pval",],3)

########################
# If using the local turn method to get at biological relevant steps
movResTA <- movLikelihoods(mov, PRdetails=TRUE, TAc=10)
movResTA

# To look at the best model compare AICc
AICcResTA <- unlist(movResTA$mle)[grep("AICc", names(unlist(movResTA$mle)))] # CRW is still the best
AICcResTA - min(AICcResTA) # By far


###
# Look at confidence intervals, the "true"/simulated value don't always fall in the 95%CI
cbind(simVal=l, movRes$CI$E,
      inInt = movRes$CI$E[,2] <= l & movRes$CI$E[,3] >= l)

# Look at the profile likelihood CI over the range from the quad approximation
# THe SL and TA are separated for models other than CCRW
layout(matrix(1:2,nrow=1))
par(mar=c(4,4,0.5,0.5), mgp=c(2.5,0.8,0))
# Step length based on BW
ciPLSL <- ciEpl(mov, movRes$mle$BW, rangePar=movRes$CI$E[,2:3])
# You'll get warnings if the range values are outside the possible range for the parameter
ciPLSL
# TA based on CRW
ciPLTA <- ciKpl(mov, movRes$mle$CRW, rangePar=movRes$CI$K[,2:3])
# You'll get warnings if the range values are outside the possible range for the parameter
ciPLTA

# Clearly the quad approximation are not perfect estimates the CI interval
rangeB <- cbind(movRes$CI$E[,2]*0.95,movRes$CI$E[,3]*1.05)
ciPLSL <- ciEpl(mov, movRes$mle$BW, rangePar=rangeB)
cbind(simVal=l, ciPLSL,
      inInt = ciPLSL[,2] <= l & ciPLSL[,3] >= l)
rangeB <- cbind(movRes$CI$K[,2]*0.95,movRes$CI$K[,3]*1.05)
ciPLTA <- ciKpl(mov, movRes$mle$CRW, rangePar=rangeB)
cbind(simVal=k, ciPLTA,
      inInt = ciPLTA[,2] <= k & ciPLTA[,3] >= k)

###
# To look at the estimates from the truncated version
# Step length based on BW
ciPLSLT <- ciTEpl(mov, movRes$mle$TBW, rangePar=movRes$CI$TE[,2:3])
ciPLSLT
# TA based on CRW
ciPLTAT <- ciKpl(mov, movRes$mle$TCRW, rangePar=movRes$CI$K[,2:3])
rbind(ciPLTAT,ciPLTA) # should be the same as long as you're looking over the same range

###
# Compare resuts to simulation values
cbind(movResTA$mle$CRW[1:3], c(l,k,a)) # Not so bad except for k
# It's expected that k would be biased because we are removing all small angles

###
# Look at test of absolute fit
# Significantly different from all model with an alpha of 0.05
# But that's likely driven by the removal of the small truning angle
round(movResTA$pseudoRes$PR["pval",],3)
# If we focuss only on the step length
round(movResTA$pseudoRes$Z["pval",],3)
# We see that all turning angle (TA) distributions are insufficient (sig. dif.)
# but the SL distribution of CCRW, E, TE are not sig. dif
