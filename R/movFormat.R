movFormat <- function(movltraj, TAc=0){
	##
	# The function input needs to be a ltraj(typeII) object
#	if (class(movltraj)[1] != "ltraj"){warning("The input needs to be of class ltraj.")}

	# ltraj object can have multiple bursts
	# but in our case we are applying the data set to individuals movement only.
	# If you have multiple burst for the same individual that you want to analyse together.
	# You have to join them before hand and setNA.

#	if(length(movltraj)>1){
#		warning("Analysis only applied to the first burst of the ltraj")
#	}

	##
	# The likelihood functions are made for locations taken at (almost) regular time intervals
	# To make sure that step length and turning angle are representing the regular time intervals
	# we use setNA() to identify missing locations
	
	# I recommend setting the appropriate time difference before using the movLikelihoods() function
	# and not use the first dt but rather use the dt that you know you want to use
	# e.g.: The grizzly bear collars are set to get location every 4hrs and they start
	# at 00:00 and I have use that information to setNA
#	movltraj <- setNA(movltraj,date.ref=movltraj[[1]]$date[1],dt=min(movltraj[[1]]$dt, na.rm=TRUE))

	# Because it is impossible to get a relative turning angle for the "steps" for which
	# the animal stayed exactly at the same location
	# (i.e. giving a value of 0 to such turning angle would bias the turning angle distribution)
	# the ltraj object remove the turning angle for consecutive
	# locations that are excatly at the same location. 
	# I have also decided to not use the 1st turning angle after the animal restart moving
	# This would have required using the step length from before the animal stop and 
	# to be conservative I've decided not to use it.
	# i.e. I've used slsp="missing" see ltraj information
#	rec(movltraj, slsp="missing")

	##
	# Although the likelihood functions are made for regular time intervals it can handle missing values.
	# Note that you need set NA for the ltraj and then remove them from the extracted SL and TA.
	# If you don't do so your relative turning angle will be wrong.
	# We only keep the steps that have both turning angle and step length
	# Because we require 3 successive locations for relative turning angle (but only 2 for step lentgh)
	# and steps at the exact smae locations have no turning angle,
	# the turning angles are the limiting.

  ########################
  # Definition based on significant turns
  # TAc is the critical turning angle in degree
  # All steps with a turning angle smaller than the crietical angle (in either direction)
  # are considered the same step
  
  # Change TAc from degree to radian
  TAc <- TAc/360*2*pi


  # The first relative turning angle should be NA
  # since it requires 3 locations to be calculated
  # But if only a subset a a movement path is taken it might not be NA
  # We set it to NA
  # There are a few subtleties here to keep in mind
  # Although you only need 2 locations for a step length,
  # you need 3 locations to make a turning angle
  movltraj[[1]]$rel.angle[1] <- NA
  
  # We are keepin only the locations associated with the step for which
  # TA are greater are equal to the critical angle or NA
  tsi <- abs(movltraj[[1]]$rel.angle) >= TAc | is.na(movltraj[[1]]$rel.angle)
  
  movltraj <- as.ltraj(movltraj[[1]][tsi,c("x","y")],
                        movltraj[[1]][tsi,"date"],
                        id(movltraj))
  
  notMisLoc <- which(!is.na(movltraj[[1]]$rel.angle))
  
  SL <- movltraj[[1]]$dist[notMisLoc]
  
  # Because I'm using both circular and CircStats
  # I need TA in both numerical format and circular format
  TA <- movltraj[[1]]$rel.angle[notMisLoc]
  TA_C <- circular(TA)
  
  # Note that the minimum and maximum values of the step lengths
  # are determine by the data.
  # So when they are used in the exponential and pareto distribution
  # they are considered as a estimated parameter.
  # See Edwards 2011 (supplementary material) for more details
  SLmin <- min(SL)
  SLmax <- max(SL)	
  
  # Need the missing locations time diff for the HMM
  missL <- diff(notMisLoc)
  
  
  notMisLoc <- notMisLoc-notMisLoc[1]+1
  
  # Length of time series
  n <- length(SL)
  
	return(list('SL'=SL,'TA_C'=TA_C, 'TA'=TA, 'SLmin'=SLmin,'SLmax'=SLmax,
		'missL'=missL,'n'=n,'notMisLoc'=notMisLoc, 'movltraj'=movltraj))
}
