library(CCRWvsLW)
l <- 0.01
k <- 10
a <- 1 
mov <- simmCRW(500,l,k,a)
movRes <- movLikelihoods(mov, PRdetails=TRUE)
movRes

# To look at the best model compare AICc
AICcRes<- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] # CRW is the best
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
AICcResTA <- unlist(movResTA$mle)[grep("AICc", names(unlist(movResTA$mle)))] # CCRW is still the best
AICcResTA - min(AICcResTA) # By far

# Compare resuts to simulation values
cbind(movResTA$mle$CRW[1:3], c(l,k,a)) # Not so bad except for k
# It's expected that k would be biased because we are removing all small angles


# Look at test of absolute fit
# Significantly different from all model with an alpha of 0.05
# But that's likely driven by the removal of the small truning angle
round(movResTA$pseudoRes$PR["pval",],3)
# If we focuss only on the step length
round(movResTA$pseudoRes$Z["pval",],3)
# We see that all turning angle (TA) distributions are insufficient (sig. dif.)
# but the SL distribution of CCRW, E, TE are not sig. dif
