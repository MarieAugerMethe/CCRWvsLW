# Simple TLW example

library(CCRWvsLW)
mu <- 2
a <- 1 
b <- 10000
mov <- simmTLW(500,mu,a,b)
movRes <- movLikelihoods(mov, PRdetails=TRUE)
movRes

# To look at the best model compare AICc
AICcRes <- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] 
AICcRes - min(AICcRes) # TLW is the best, followed by LW (which should be considered more or less the same model)

# Compare resuts to simulation values
cbind(movRes$mle$TLW[1:3], c(mu,a,b)) 
# Pretty good, 
# except for b, 
# but that's because you need a really big sample size to have a chane to have really long step

###
# Look at the profile likelihood CI over the range from the quad approximation
par(mar=c(4,4,0.5,0.5), mgp=c(2.5,0.8,0))
ciPL <- ciTLWpl(mov, movRes$mle$TLW, rangePar=movRes$CI$TLW[,2:3])
# You'll get warnings if the range values are outside the possible range for the parameter
ciPL
# Clearly the quad approximation are not perfect estimates the CI interval
rangeB <- cbind(movRes$CI$TLW[,2]*0.95,movRes$CI$TLW[,3]*1.05)
ciPL <- ciTLWpl(mov, movRes$mle$TLW, rangePar=rangeB)
cbind(simVal=mu, ciPL,
      inInt = ciPL[,2] <= mu & ciPL[,3] >= mu)
# Still simulated value not always in CI

###
# Look at test of absolute fit
# With an alpha of 0.05, not significantly different from LW & TLW
round(movRes$pseudoRes$PR["pval",],3)

########################
# If using the local turn method to get at biological relevant steps
movResTA <- movLikelihoods(mov, PRdetails=TRUE, TAc=10)
movResTA

# To look at the best model compare AICc
AICcResTA <- unlist(movResTA$mle)[grep("AICc", names(unlist(movResTA$mle)))]
AICcResTA - min(AICcResTA) # TLW and LW are stil the best

# Compare resuts to simulation values
cbind(movResTA$mle$TLW[1:3], c(mu,a,b))  # Similar results

# Look at test of absolute fit
# Significantly different from all model with an alpha of 0.05
# But that's likely driven by the removal of the small turning angle
round(movResTA$pseudoRes$PR["pval",],3)
# If we focuss only on the step length
round(movResTA$pseudoRes$Z["pval",],3)
# We see that all turning angle (TA) distributions are insufficient (sig. dif.)
# but the SL distribution of LW & TLW are not sig. dif
