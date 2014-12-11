# Call the EM algorithm function written in Rcpp and set default values
emHMM <- function(SL, TA, missL, SLmin, lambda, gamm, delta=c(0.5,0.5), kapp, notMisLoc,
                        maxiter=10000, tol=1e-5){
  res <- EMHMM(SL, TA, missL, notMisLoc, SLmin, lambda, gamm, delta, kapp, maxiter,
                tol, optK)
  return(res)
}

# To function needed to numerically get the MLE of kappa in the EM algorithm
optK <- function(t2){
  mk <- optimize(function(k) abs((besselI(k,1)/besselI(k,0))-t2),
                 c(0.0000001,709))$minimum
  return(mk)
} 


HMM.lalphabeta <- function(SL,TA,missL,SLmin=min(SL),lambda, kapp, gamma, delta, notMisLoc){

	# This is the real length of the time series (include missing data)
	n <- length(SL) + sum(missL-1) # or maybe just notMisLoc[length(notMisLoc)]

	###
	# Initialising the alpha and beta vectors
	# and creating a matrix will holds alpha and beta values for all t
	# Note that to limit the underflow problem we are evaluating them as log values 
	lalpha <- lbeta <- matrix(NA,2,n)

	# Making a matrix with all the probability of observations for each state
	# Note that for missing location observation probablility is 1
	obsProb <- matrix(1,n,2)
	obsProb[notMisLoc,1] <- dexp((SL-SLmin),lambda[1]) * dvm(TA, 0, 0)
	obsProb[notMisLoc,2] <- dexp((SL-SLmin),lambda[2]) * dvm(TA, 0, kapp)

	# According to Zucchini and MacDonald (2009):
	# alpha_t = alpha_{t-1} %*% gamma %*%P(x_t)
	# alpha_1 = delta %*% P(x_1)
	# phi_t = alpha_t / w_t
	# w_t = sum(alpha_t)
	foo <- delta*obsProb[1,] # delta_0 %*% P(x_1) = alpha_1
	sumfoo <- sum(foo) # w_1
	lscale <- log(sumfoo) # log(w_1)
	foo <- foo/sumfoo # alpha_1/w_1 = phi_1
	lalpha[,1] <- log(foo) + lscale # log(phi_1) + log(w_1) = log(alpha_1)
		
	# starting the log alpha recursion
	for (i in 2:n){
		foo <- foo%*% gamma *obsProb[i,] # phi_{i-1}%*%gamma%*%P(x_{i}) = alpha_i/w_{i-1}
		sumfoo <- sum(foo) # sum(alpha_i/w_{i-1}) = w_i/w_{i-1}
		lscale <- lscale+log(sumfoo) #log(w_{i-1})+log(w_i)-log(w_{i-1})=log(w_i)
		foo <- foo/sumfoo # (alpha_i/w_{i-1})/(w_i/w_{i-1}) = alpha_i/w_{i} =phi_{i}
		lalpha[,i] <- log(foo)+lscale # log(phi_i) + log(w_i) = log(alpha_i)
	}

	# we also want the backward probability to get the expected values
	# tr(beta_t) = gamma %*% obsProb_{t+1} %*% tr(beta{t+1}) %*% tr(1)
	# I think similarly to alpha the betas are weighted
	# psi_t = beta_t/w_t
	# w_t = beta_t %*% tr(1)
	lbeta[,n] <- rep(0,2) # beta_T = 1 so log(beta_T)=0 
	foo <- rep(0.5,2) # beta_T/w_T =psi_T 
	lscale <- log(2) # log(w_T)
	for (i in (n-1):1){
		foo <- gamma %*% (obsProb[i+1,]*foo) # gamma%*%P(x_{t+1}))*psi_{t+1} = beta_t/w_{t+1}
		lbeta[,i] <-log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
		sumfoo <- sum(foo) # w_t/w_{t+1}
		foo <- foo/sumfoo # (beta_t/w_{t+1})/(w_t/w_{t+1}) = psi_t
		lscale <- lscale + log(sumfoo) # log(w_{t+1}) + log(w_t) - log(w_{t+1}) = log(w_t)
	}
	list(la=lalpha, lb=lbeta)
}

# To be used in the pseudo-residual function
# (and for a graph that show the probability of being in each behavior)
# The weight function kind of calculates the probability of being in either behaviour
# based on all observations other than the one at time t
# It is define as: w_i(t) = d_t(t)/(sum_{j=1}^{m}{dj(t)})
# where d_i(t) = (alpha_{t-1}%*%gamma)[i] * beta_t(i)
# The weights function and weights are using for both SL and TA
# and uses both SL and TA
HMMwi <- function(SL,TA,missL,SLmin,lambda,kapp,gamma,delta,notMisLoc){
  # This calculates the weights for the pseudo-residuals
  n <- length(SL) + sum(missL-1)
  fb <- HMM.lalphabeta(SL,TA,missL,SLmin,lambda, kapp, gamma, delta, notMisLoc)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)
  w <- matrix(NA,2,n)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%gamma)*
      exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  return(w)
}