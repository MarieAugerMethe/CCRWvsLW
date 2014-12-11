#################################
# Created by: Marie Auger-Methe
# Date: May 24, 2011
# Updated: April 8, 2013

#######################################
# CCRW_HMM
simmCCRW <- function(n, gII, gEE, lI, lE, kE, a, dI=NULL){
  n <- n+1 # Because the ltraj object places NA at the first and last location
  
  lambda <- c(lI,lE)
  ai2 <- c((1-gII),gEE)
  
  if(is.null(dI)){
    #######################
    # Calculate the stationary distribution of the Markov Chain
    tpm <- matrix(c(gII, (1-gEE), (1-gII), gEE), nrow=2, byrow=TRUE)
    dI <- eigen(tpm)$vectors[,1]
    # Need to normalise the eigen vector to sum to 1
    dI <- dI/sum(dI)
    dI <- dI[1]
  }
  
  pc <- sslta(n, a, lambda, kE, ai2, (1-dI))
  
  coord <- p2c(pc[,1],pc[,2])
  idName <- paste("CCRWgII", as.character(round(gII,3)),
                  "gEE", as.character(round(gEE,3)),sep="")
  dat <- 1:(n+1)
  class(dat) <- c("POSIX", "POSIXct")
  mov <- as.ltraj(coord, dat, id=idName)
  
  return(mov)
}


#######################################
# LW

simmLW <- function(n, mu, a){
# This functions simmulated a Levy walk (LW)
# Using the pareto distribution (power law) for the step length distribution
# and a von Mises distribution for the the turning angle distribution.
# The von Mises distribution has the kappa set to 0 which reduces it
# to an uniform distribution.
# The turning anle is the relative turning angle.
# The function returns a ltraj object (from adehabitat package)

  n <- n+1 # Because the ltraj object places NA at the first and last location
  
	###########################
	# Parameters

	# Specified as input
	# n is the number of steps
	# mu is the shape parameter of the Pareto
	# a is location parameter (the start of the tail of the power law)
	# Note that mu is the equivalent of c+1 in Forbes et al 2011 description of 
	# the Pareto distribution

	# Defined in the LW
	# Uniform turning angle distribution: mu =0, kappa=0

	##########################
	# Sampling from the step length and turn angle distributions

	# The power law distribution used in the Levy walk
	# is the pareto distribution
	# However the way the parameters are written in the animal movement litterature
	# differs from the way the Pareto distribution is defined in R
	# The difference is:
	# k <- mu-1
	# To test:
	# dpareto(stepL, a, (mu-1))
	# ((mu-1)/(a^(1-mu))*stepL^(-mu))

	SL <- rpareto(n, a, (mu-1))
	TA <- runif(n, 0, 2*pi)
	
	########################
	# Creating ltraj object
	
  # Change the step length and relative turn angle 
  # into cartesian coordinates
	coord <- p2c(SL,TA)
  
	# Dates & time interval
	dat <- 1:(n+1)
	class(dat) <- c("POSIX", "POSIXct")
  
	idName <- paste("LWmu", as.character(round(mu,3)), sep="")
  
	mov <- as.ltraj(coord, dat, id=idName)
	return(mov)
}

#######################################
# TLW

simmTLW <- function(n, mu, a, b){

  n <- n+1 # Because the ltraj object places NA at the first and last location
  
# This functions simmulated a truncated Levy walk (TLW)
# Using a truncated pareto distribution (truncated power law)
# for the step length distribution
# and a von Mises distribution for the the turning angle distribution.
# The von Mises distribution has the kappa set to 0 which reduces it
# to an uniform distribution.
# The turning angle is the relative turning angle.
# The function returns a ltraj object (from adehabitat package)


	###########################
	# Parameters

	# Specified as input
	# n is the number of steps
	# mu is the shape parameter of the truncated power law 
	# a is the location parameter (the start of the tail of the power law)
	# b is the upper bound of the power law

	# Defined in the TLW model
	# Uniform turning angle distribution: mu=0, kappa=0
	
	##########################
	# Sampling from the step length and turn angle distributions

	# The truncated power law distribution used for the TLW
	# is the truncated pareto distribution
	# However the way the parameters are written
	# in the animal movement litterature
	# differs from the way the Pareto distribution is defined in R
	# The difference is: 
	# k <- mu-1
	# To test:
	# dtpareto(stepL, a, b, (mu-1))
	# ((mu-1)/(a^(1-mu)-b^(1-mu))*stepL^(-mu))

	SL <- rtruncpareto(n, a, b, (mu-1)) # Used to be rtpareto
	TA <- runif(n, 0, 2*pi)
	
  ######################
  # Creating an ltraj object
	# Change the step length and relative turn angle 
	# into cartesian coordinates
	coord <- p2c(SL,TA)
	
	# Dates & time interval
	dat <- 1:(n+1)
	class(dat) <- c("POSIX", "POSIXct")
	
	idName <- paste("TLWmu", as.character(round(mu,3)), sep="")
	mov <- as.ltraj(coord, dat, id=idName)
	
	return(mov)
}

#######################################
# BW
simmBW <- function(n, l, a){
  # This functions simmulated a Brownian walk
  # Using a exponential distribution for the step length distribution
  # and a von Mises distribution for the the turning angle distribution
  # The step length smaple from the exponential distribution
  # is added to the minimum step length, a. 
  # The von Mises distribution has the kappa set to 0 which reduces it
  # to an uniform distribution.
  # The turning anle if the relative turning angle.
  # The function returns a ltraj object (from adehabitat package)
  n <- n+1 # Because the ltraj object places NA at the first and last location
  
  ###########################
  # Parameters
  
  ##
  # Specified as input
  # n is the number of steps
  # lambda is the rate parameter (the inverse of the scale paramter) of the exponetial distribution
  # a is the minimum step length
  
  # Defined in the BW
  # Uniform turning angle: mu =0, kappa =0
  
  ##########################
  # Sampling from the step length and turn angle distributions
  SL <- rexp(n, l) + a
  TA <- runif(n, 0, 2*pi)
  
  
  ######################
  # Creating an ltraj object
  # Change the step length and relative turn angle 
  # into cartesian coordinates
  coord <- p2c(SL,TA)
  
  # Dates & time interval
  dat <- 1:(n+1)
  class(dat) <- c("POSIX", "POSIXct")
  
  idName <- paste("BWl", as.character(round(l,3)), sep="")
  mov <- as.ltraj(coord, dat, id=idName)
  
  return(mov)
}

#######################################
# CRW

simmCRW <- function(n, l, k, a){
# This function simmulates a correlated random walk (CRW)
# Using an exponential distribution for the step length distribution
# and a von Mises distribution for the the turning angle distribution.
# The step length sample from the esponential distribution are added to
# the minimum step length.
# The von Mises distribution is centered around 0 to represent foward persistance
# The turning angle is the relative turning angle.
# The function returns a ltraj object (from adehabitatLT package)
  n <- n+1 # Because the ltraj object places NA at the first and last location

	###########################
	# Parameters

	##
	# Specified as input
	# n is the number of steps
	# lambda is the rate parameter (the inverse of the scale paramter) of the exponetial distribution
	# kapp is the concentration parameter for the turning angle representing how "correlated" the random walk is
	# a is the minimum step length

	##
	# Defined in the correlated random walk model
	# To represent the forward persistence mu is set to 0

	##########################
	# Sampling from the step length and turn angle distributions
	SL <- rexp(n, l)+a
	TA <- rvm(n, 0, k)

	######################
	# Creating an ltraj object
	# Change the step length and relative turn angle 
	# into cartesian coordinates
	coord <- p2c(SL, TA)
	
	# Dates & time interval
	dat <- 1:(n+1)
	class(dat) <- c("POSIX", "POSIXct")
	
	idName <- paste("CRWl", as.character(round(l,3)),
                  "k", as.character(round(k,3)), sep="")
	mov <- as.ltraj(coord, dat, id=idName)

	return(mov)
}
