# Example of simulating and fitting a CCRW

library(CCRWvsLW)
gII <- 0.9
gEE <- 0.9
lI <- 0.01
lE <- 0.001
kE <- 10
a <- 1 
mov <- simmCCRW(500,gII,gEE,lI,lE,kE,a,0.5)
movRes <- movLikelihoods(mov, PRdetails=TRUE)
movRes

# To look at the best model compare AICc
AICcRes<- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] # CCRW is the best
AICcRes - min(AICcRes) # By far

# Compare parameter estimates to simulation values
cbind(simVal=c(gII,gEE,lI,lE,kE), parEst=movRes$mle$CCRW[1:5]) # Pretty good

# Look at confidence intervals, the "true"/simulated value don't always fall in the 95%CI
cbind(simVal=c(gII,gEE,lI,lE,kE), movRes$CI$CCRW,
      inInt = movRes$CI$CCRW[,2]<=c(gII,gEE,lI,lE,kE) & movRes$CI$CCRW[,3]>=c(gII,gEE,lI,lE,kE))

# Look at the profile likelihood CI over the range from the quad approximation
par(mar=c(4,4,0.5,0.5), mgp=c(2.5,0.8,0))
ciPL <- ciCCRWpl(mov, movRes$mle$CCRW, rangePar=movRes$CI$CCRW[,2:3])
# You'll get warnings if the range values are outside the possible range for the parameter
ciPL
# Clearly the quad approximation are underestimating the CI interval
rangeB <- cbind(movRes$CI$CCRW[,2]*0.95,movRes$CI$CCRW[,3]*1.05)
ciPL <- ciCCRWpl(mov, movRes$mle$CCRW, rangePar=rangeB)
cbind(simVal=c(gII,gEE,lI,lE,kE), ciPL,
      inInt = ciPL[,2]<=c(gII,gEE,lI,lE,kE) & ciPL[,3]>=c(gII,gEE,lI,lE,kE))
# Now simulated value is in the CI

# Look at test of absolute fit
# With an alpha of 0.05, not significantly different from CCRW
round(movRes$pseudoRes$PR["pval",],3)

########################
# If using the local turn method to get at biological relevant steps
movResTA <- movLikelihoods(mov, PRdetails=TRUE, TAc=10)
movResTA

# To look at the best model compare AICc
AICcResTA <- unlist(movResTA$mle)[grep("AICc", names(unlist(movResTA$mle)))] # CCRW is still the best
AICcResTA - min(AICcResTA) # By far

# The difference with the threshold angle is smaller (likelily in part because the data set is smaller), 
# but the delta AIC is still huge
cbind(AICcRes,AICcResTA)

# Compare resuts to simulation values
cbind(movResTA$mle$CCRW[1:5], c(gII,gEE,lI,lE,kE)) # Not so bad except for kE
# It's expected that kE would be biased because we are removing all small angles

# Look at test of absolute fit
# Significantly different from all model with an alpha of 0.05
# But that's likely driven by the removal of the small truning angle
round(movResTA$pseudoRes$PR["pval",],3)
# If we focuss only on the step length
round(movResTA$pseudoRes$Z["pval",],3)
# We see that all turning angle (TA) distributions are insufficient (sig. dif.)
# but the SL distribution of CCRW is not sig. dif
