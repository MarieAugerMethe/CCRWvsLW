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
	n <- length(SL) + sum(missL-1)

	###
	# Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
	# Note that to limit the underflow problem we are evaluating them as log values.
	lalpha <- lbeta <- matrix(NA,2,n)

	# Making a matrix with all the probability of observations for each state.
	# Note that for missing location observation probablility is 1.
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

## for ww
HMMww.lalphabeta <- function(SL,TA, missL, sh, sc, rE, gamma, delta, notMisLoc){
  
  # This is the real length of the time series (include missing data)
  n <- length(SL) + sum(missL-1)
  
  ###
  # Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
  # Note that to limit the underflow problem we are evaluating them as log values. 
  lalpha <- lbeta <- matrix(NA,2,n)
  
  # Making a matrix with all the probability of observations for each state
  # Note that for missing location observation probablility is 1
  obsProb <- matrix(1,n,2)
  obsProb[notMisLoc,1] <-  dweibull(SL, sh[1], sc[1]) * dwrpcauchy(TA, 0, 0)
  obsProb[notMisLoc,2] <-  dweibull(SL, sh[2], sc[2]) * dwrpcauchy(TA, 0, rE)
  
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
    foo <- foo%*% gamma * obsProb[i,] # phi_{i-1}%*%gamma%*%P(x_{i}) = alpha_i/w_{i-1}
    sumfoo <- sum(foo) # sum(alpha_i/w_{i-1}) = w_i/w_{i-1}
    lscale <- lscale + log(sumfoo) #log(w_{i-1})+log(w_i)-log(w_{i-1})=log(w_i)
    foo <- foo/sumfoo # (alpha_i/w_{i-1})/(w_i/w_{i-1}) = alpha_i/w_{i} =phi_{i}
    lalpha[,i] <- log(foo) + lscale # log(phi_i) + log(w_i) = log(alpha_i)
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
    lbeta[,i] <- log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
    sumfoo <- sum(foo) # w_t/w_{t+1}
    foo <- foo/sumfoo # (beta_t/w_{t+1})/(w_t/w_{t+1}) = psi_t
    lscale <- lscale + log(sumfoo) # log(w_{t+1}) + log(w_t) - log(w_{t+1}) = log(w_t)
  }
  list(la=lalpha, lb=lbeta)
}

## for HSMM
HSMM.lalphabeta <- function(SL,TA, missL, sh, sc, rE, gamSize, gamPr, notMisLoc,m = c(10,10)){
  
  gamma <- gen.Gamma.repar(m,gamSize,gamPr) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  
  # This is the real length of the time series (include missing data)
  n <- length(SL) + sum(missL-1)
  
  obsProb <- matrix(rep(1,sum(m)*n),nrow=n)
  # Making a matrix with all the probability of observations for each state
  # Note that for missing location observation probablility is 1
  obsProb[notMisLoc,1:m[1]] <-  dweibull(SL, sh[1], sc[1]) * dwrpcauchy(TA, 0, 0)
  obsProb[notMisLoc,(m[1]+1):sum(m)] <-  dweibull(SL, sh[2], sc[2]) * dwrpcauchy(TA, 0, rE)
  
  ###
  # Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
  # Note that to limit the underflow problem we are evaluating them as log values. 
  lalpha <- lbeta <- matrix(NA,sum(m),n)
  
  
  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*obsProb[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  
  
  # we also want the backward probability to get the expected values
  # tr(beta_t) = gamma %*% obsProb_{t+1} %*% tr(beta{t+1}) %*% tr(1)
  # I think similarly to alpha the betas are weighted
  # psi_t = beta_t/w_t
  # w_t = beta_t %*% tr(1)
  lbeta[,n] <- rep(0,sum(m)) # beta_T = 1 so log(beta_T)=0 
  foo <- rep(0.5,sum(m)) # beta_T/w_T =psi_T 
  lscale <- log(sum(m)) # log(w_T)
  for (i in (n-1):1){
    foo <- gamma %*% (obsProb[i+1,]*foo) # gamma%*%P(x_{t+1}))*psi_{t+1} = beta_t/w_{t+1}
    lbeta[,i] <- log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
    sumfoo <- sum(foo) # w_t/w_{t+1}
    foo <- foo/sumfoo # (beta_t/w_{t+1})/(w_t/w_{t+1}) = psi_t
    lscale <- lscale + log(sumfoo) # log(w_{t+1}) + log(w_t) - log(w_{t+1}) = log(w_t)
  }
  list(la=lalpha, lb=lbeta)
}

## for HSMM
HSMMp.lalphabeta <- function(SL,TA, missL, sh, sc, rE, gamma, notMisLoc,m = c(10,10)){
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  
  # This is the real length of the time series (include missing data)
  n <- length(SL) + sum(missL-1)
  
  obsProb <- matrix(rep(1,sum(m)*n),nrow=n)
  # Making a matrix with all the probability of observations for each state
  # Note that for missing location observation probablility is 1
  obsProb[notMisLoc,1:m[1]] <-  dweibull(SL, sh[1], sc[1]) * dwrpcauchy(TA, 0, 0)
  obsProb[notMisLoc,(m[1]+1):sum(m)] <-  dweibull(SL, sh[2], sc[2]) * dwrpcauchy(TA, 0, rE)
  
  ###
  # Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
  # Note that to limit the underflow problem we are evaluating them as log values. 
  lalpha <- lbeta <- matrix(NA,sum(m),n)
  
  
  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*obsProb[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  
  
  # we also want the backward probability to get the expected values
  # tr(beta_t) = gamma %*% obsProb_{t+1} %*% tr(beta{t+1}) %*% tr(1)
  # I think similarly to alpha the betas are weighted
  # psi_t = beta_t/w_t
  # w_t = beta_t %*% tr(1)
  lbeta[,n] <- rep(0,sum(m)) # beta_T = 1 so log(beta_T)=0 
  foo <- rep(0.5,sum(m)) # beta_T/w_T =psi_T 
  lscale <- log(sum(m)) # log(w_T)
  for (i in (n-1):1){
    foo <- gamma %*% (obsProb[i+1,]*foo) # gamma%*%P(x_{t+1}))*psi_{t+1} = beta_t/w_{t+1}
    lbeta[,i] <- log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
    sumfoo <- sum(foo) # w_t/w_{t+1}
    foo <- foo/sumfoo # (beta_t/w_{t+1})/(w_t/w_{t+1}) = psi_t
    lscale <- lscale + log(sumfoo) # log(w_{t+1}) + log(w_t) - log(w_{t+1}) = log(w_t)
  }
  list(la=lalpha, lb=lbeta)
}

## for HSMM with von Mises and exponential
HSMMpo.lalphabeta <- function(SL,TA, missL, lam, kE, gamma, notMisLoc,m = c(10,10)){
  SLmin <- min(SL)
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  
  # This is the real length of the time series (include missing data)
  n <- length(SL) + sum(missL-1)
  
  obsProb <- matrix(rep(1,sum(m)*n),nrow=n)
  # Making a matrix with all the probability of observations for each state
  # Note that for missing location observation probablility is 1
  obsProb[notMisLoc,1:m[1]] <-  dexp(SL-SLmin, lam[1]) * dvm(TA, 0, 0)
  obsProb[notMisLoc,(m[1]+1):sum(m)] <-  dexp(SL-SLmin, lam[2]) * dvm(TA, 0, kE)
  
  ###
  # Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
  # Note that to limit the underflow problem we are evaluating them as log values. 
  lalpha <- lbeta <- matrix(NA,sum(m),n)
  
  
  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*obsProb[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  
  
  # we also want the backward probability to get the expected values
  # tr(beta_t) = gamma %*% obsProb_{t+1} %*% tr(beta{t+1}) %*% tr(1)
  # I think similarly to alpha the betas are weighted
  # psi_t = beta_t/w_t
  # w_t = beta_t %*% tr(1)
  lbeta[,n] <- rep(0,sum(m)) # beta_T = 1 so log(beta_T)=0 
  foo <- rep(0.5,sum(m)) # beta_T/w_T =psi_T 
  lscale <- log(sum(m)) # log(w_T)
  for (i in (n-1):1){
    foo <- gamma %*% (obsProb[i+1,]*foo) # gamma%*%P(x_{t+1}))*psi_{t+1} = beta_t/w_{t+1}
    lbeta[,i] <- log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
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
  fb <- HMM.lalphabeta(SL, TA, missL, SLmin, lambda, kapp, gamma, delta, notMisLoc)
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

# To be used in the pseudo-residual function - ww
# (and for a graph that show the probability of being in each behavior)
HMMwiww <- function(SL, TA, missL, sh, sc, rE, gamma, delta, notMisLoc){
  # This calculates the weights for the pseudo-residuals
  n <- length(SL) + sum(missL-1)
  fb <- HMMww.lalphabeta(SL, TA, missL, sh, sc, rE, gamma, delta, notMisLoc)
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



# for hsm
HSMMwi<- function(SL, TA, missL, notMisLoc, gamSize, gamPr, sc, sh, rE, m=c(10,10)){
  gamma <- gen.Gamma.repar(m,gamSize,gamPr) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  # This calculates the weights for the pseudo-residuals
  n <- length(SL) + sum(missL-1)
  fb <- HSMM.lalphabeta(SL, TA, missL, sh, sc, rE, gamSize, gamPr, notMisLoc)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)
  w <- matrix(NA,sum(m),n)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%gamma)*
      exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  w <- rbind(colSums(w[1:m[1],]),colSums(w[(m[1]+1):sum(m),]))
  return(w)
}


# for hsm with poisson
HSMMpwi<- function(SL, TA, missL, notMisLoc, lamb, sc, sh, rE, m=c(10,10)){
  gamma <- gen.Gamma.pois(m,lamb) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  # This calculates the weights for the pseudo-residuals
  n <- length(SL) + sum(missL-1)
  fb <- HSMMp.lalphabeta(SL, TA, missL, sh, sc, rE, gamma, notMisLoc)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)
  w <- matrix(NA,sum(m),n)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%gamma)*
      exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  w <- rbind(colSums(w[1:m[1],]),colSums(w[(m[1]+1):sum(m),]))
  return(w)
}

# for hsm with poisson and von Mises and exponential
HSMMpowi<- function(SL, TA, missL, notMisLoc, lamb, lambda, kE, m=c(10,10)){
  gamma <- gen.Gamma.pois(m,lamb) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  # This calculates the weights for the pseudo-residuals
  n <- length(SL) + sum(missL-1)
  fb <- HSMMpo.lalphabeta(SL, TA, missL, lambda, kE, gamma, notMisLoc)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)
  w <- matrix(NA,sum(m),n)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%gamma)*
      exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  w <- rbind(colSums(w[1:m[1],]),colSums(w[(m[1]+1):sum(m),]))
  return(w)
}

# This is the actual EM algorithm but formatted for the confidence interval
EM_CCRW_HMM_CI <- function(SL,TA_N,x,missL,notMisLoc,parF,maxiter=300,tol=5e-5){
  n <- length(SL) + sum(missL-1)
  # Note I have made the maxiter smaller for CI
  gII <- x[1]
  gEE <- x[2]
  lI <- x[3]
  lE <- x[4]
  kE <- x[5]
  dI <- x[6]
  dE <- 1-x[6]
  
  ##
  # To allow to fix some of the parameters
  # This is need for the profile likelihood CI
  if(length(parF)>0){
    for (i in 1:length(parF)){
      assign(names(parF[i]),parF[[i]])
    }
  }
  
  gamma <- matrix(c(gII,1-gEE,1-gII,gEE),nrow=2)
  lambda <- c(lI,lE)
  kapp <- kE
  delta <- c(dI,dE)
  
  # To initialise the EM algorithm you need to choose a set of parameters
  lambda.next <- lambda
  gamma.next <- gamma
  delta.next <- delta
  kappa.next <- kapp
  
  
  # Making a matrix that will store all the log probability of observations for each state
  lobsProb <- matrix(0,nrow=n,ncol=2)
  # Function that calculates the log probabilities
  lop <- function(){
    # State 1: F
    # dexp((SL-SLmin)|lambda_F) * dvonmises(TA|mu_F,kappa_F)
    lobsProb[notMisLoc,1] <-   dexp((SL-SLmin),lambda[1],log=TRUE) + log(dvm(TA_N, 0, 0))
    # State 2: T
    # dexp((SL-SLmin)|lambda_T) * dvonmises(TA|mu_T,kappa_T)
    lobsProb[notMisLoc,2] <-   dexp((SL-SLmin),lambda[2],log=TRUE) + log(dvm(TA_N, 0, kapp))
    return(lobsProb)
  }
  
  for (iter in 1:maxiter){
    lprobObs <- lop() # probability of all observations
    fb <- HMM.lalphabeta(SL,TA_N,missL,SLmin,lambda, kapp, gamma, delta, notMisLoc) # foward and backward prob
    la <- fb$la
    lb <- fb$lb
    # Use alpha to get the likelihood
    # remember that LT = alpha_T %*%tr(1)
    # Using c is just to limit underflow
    c <- max(la[,n])
    
    if(is.nan(c)==T){
      warning(paste("Underflow problem cannot calculate the log likelihood.",
                    "This is iteration", iter))
      return(list(lambda=NA, gamma=NA, delta=NA, kapp=NA, mllk=NA))
    }
    
    llk <- c + log(sum(exp(la[,n]-c))) # log(max(alpha_T)) + log(sum(alpha_T/max(alpha_T))) 
    
    # Maximum likelihood estimate based on complete data log likelihood (CDLL)
    for(j in 1:2){
      for(k in 1:2){
        
        # Calculating: sum_{t=2}^T(E(v_{jk}(t)))
        # which is needed for the mle of gamma (see below)
        # Expected binary transition value
        # E(v_jk(t))=a_{t-1}(j)*gamma_{jk}*p_k(x_t)*b_t(k)/L_T
        
        gamma.next[j,k] <- gamma[j,k] * # gamma_{jk}
          sum(exp(
            la[j,1:(n-1)] +  # la_{1:(T-1)}(j))
              lprobObs[2:n,k] + # lp_k(2:T)
              lb[k,2:n] - # lb_k(2:T)
              llk)) 
      }
      # MLE of CDLL for lambda
      # sum_{t=1}^T(E(u_j(t)) /sum_{t=1}^T(E(u_j(t)*SL_t)
      # Expected binary state value
      # E(u_j(t))=a_t(j)*b_t(j)/L_T
      lambda.next[j] <- sum(
        exp(la[j,notMisLoc]+lb[j,notMisLoc]-llk) # E(u_j(t))
      )/
        sum(exp(la[j,notMisLoc]+lb[j,notMisLoc]-llk)*(SL-SLmin))			
      
    }
    # MLE of CDLL for kappa
    # I_1(kappa)/I_0(kappa) = sum (E(u_j(t))*cos(x_t))/sum(E(u_j(t)))
    # Where I_1 and I_0 are the modified bessel functions
    # To get kappa
    # 0 = I_1(kappa)/I_0(kappa) - sum (E(u_j(t))*cos(x_t))/sum(E(u_j(t)))
    # It's is possible that this equation is biased for kappa
    # Although my simulation should verify this potential problem
    # I can use optimize to get the kappa that minimize the function above
    # For foraging behaviour the kappa is assumed to be 0 (uniform distribution)
    # so only for travelling behaviour So j=2
    
    term2 <- sum(exp(la[2,notMisLoc]+lb[2,notMisLoc]-llk)*cos(TA_N))/
      sum(exp(la[2,notMisLoc]+lb[2,notMisLoc]-llk))
    toMin <- function(k){
      abs((besselI(k,1)/besselI(k,0))-term2)
    }
    minKappa <- optimize(toMin,c(0.0000001,709)) # if Kappa greater than 709 besselI(k,1) = Inf
    kappa.next <- minKappa$minimum
    
    # gamma_{jk} = sum_{t=2}^T(E(v_{jk}(t)))/sum_{k=1}^m(sum_{t=2}^T(E(v_{jk}(t))))
    gamma.next <- gamma.next/rowSums(gamma.next) # apply(,1,sum) : sum row
    
    # MLE for delta
    # E(u_j(1))=a_t(1)*b_t(1)/L_T
    # delta_j = E(u_j(1))/sum_{j=1}^1(E(u_j(1))) = E(u_j(1))
    delta.next <- exp(la[,1]+lb[,1]-llk)
    delta.next <- delta.next/sum(delta.next)
    
    
    a11 <- gamma.next[1]
    a22 <- gamma.next[4]
    lambda_F <- lambda.next[1]
    lambda_T <- lambda.next[2]
    kappa_T <- kappa.next
    delta_F <- delta.next[1]
    delta_T <- delta.next[2]
    
    ##
    # To allow to fix some of the parameters
    # This is need for the profile likelihood CI
    if(length(parF)>0){
      for (i in 1:length(parF)){
        assign(names(parF[i]),parF[[i]])
      }
    }
    
    gamma.next <- matrix(c(gII,1-gEE,1-gII,gEE),nrow=2)
    lambda.next <- c(lI,lE)
    kappa.next <- kE
    delta.next <- c(dI,dE)
    
    crit <- sum(abs(lambda -lambda.next)) + # for convergence
      sum(abs(gamma -gamma.next)) +
      sum(abs(delta -delta.next)) +
      sum(abs(kapp - kappa.next))
    
    if(crit <tol){
      return(list(lambda=lambda, gamma=gamma, delta=delta, kapp=kapp, mllk=-llk))
    }
    
    # store the values to be compare on the next iteration of the loop
    
    lambda <- lambda.next
    gamma <- gamma.next
    delta <- delta.next
    kapp <- kappa.next
    
    
    # I think there are problems of unbounded likelihood
    # See Zucchini and MacDonald (2009) p.10 & 50
    # My quick fix is to stop when the estimated lambdas are Inf
    if(length(which(lambda==Inf))>0){
      warning(paste("Potential unbounded likelihood problem.", 
                    "Lambda_I = ",  lambda[1],
                    "Lambda_E = ", signif(lambda[2],4)))
      return(list(lambda=NA, gamma=NA, delta=NA, kapp=NA, mllk=NA))
    }
  }
  warning(paste("No convergence after", maxiter,"iterations"))
  return(list(lambda=NA, gamma=NA, delta=NA, kapp=NA, mllk=NA))
  
}
