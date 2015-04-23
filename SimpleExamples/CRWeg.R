# Example of simulating and fitting a CRW

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
