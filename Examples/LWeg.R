library(CCRWvsLW)
mu <- 2
a <- 1 

mov <- simmLW(500,mu,a)
plot(mov)
movRes <- movLikelihoods(mov, PRdetails=TRUE)
# You'll most likeliy have warnings, 
# becasue the HMM will have a hard time fitting a model with extremely long step length
# In this case it make, sense since the CCRW/HMM is not a good model for this LW dataset

# To look at the best model compare AICc
AICcRes <- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] 
AICcRes - min(AICcRes) # LW & TLW are the bests

# Compare resuts to simulation values
cbind(movRes$mle$TLW[1:2], c(mu,a)) # Pretty good, 

# Look at test of absolute fit
# With an alpha of 0.05, not significantly different from LW and TLW
round(movRes$pseudoRes$PR["pval",],3)

########################
# If using the local turn method to get at biological relevant steps
movResTA <- movLikelihoods(mov, PRdetails=TRUE, TAc=10)
# Again, high probability of warnings

# To look at the best model compare AICc
AICcResTA <- unlist(movResTA$mle)[grep("AICc", names(unlist(movResTA$mle)))]
AICcResTA - min(AICcResTA) # TLW and LW are stil the best

# Compare resuts to simulation values
cbind(movResTA$mle$LW[1:2], c(mu,a))  # Similar results


# Look at test of absolute fit
# Significantly different from all model with an alpha of 0.05
# But that's likely driven by the removal of the small truning angle
round(movResTA$pseudoRes$PR["pval",],3)
# If we focuss only on the step length
round(movResTA$pseudoRes$Z["pval",],3)
# We see that all turning angle (TA) distributions are insufficient (sig. dif.)
# but the SL distribution ofLW & TLW are not sig. dif
