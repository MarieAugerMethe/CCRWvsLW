######################
# Example to fit the models to a Grizzly movement track
# This is to reproduce the analysis of grizzle G2 in 
# Auger-Methe et al. (2016) Evaluating random search strategies in three mammals
# from distinct feeding guils. Journal of Animal Ecology

# The results can be compared to those of the supporting information

library(CCRWvsLW)

##
# Load data
grizzly <- read.csv("SampleData/GrizzlySample.csv", stringsAsFactors = FALSE)


##
# Transform data into ltraj object
# We have the movement in the x and y direction
x <- c(0,cumsum(grizzly$dx)[-length(grizzly$dx)])
y <- c(0,cumsum(grizzly$dy)[-length(grizzly$dy)])
grizzlyDate <- as.POSIXct(grizzly$Date,tz="UTC")

# Make ltraj object
cMov <- as.ltraj(xy=cbind(x,y), date=grizzlyDate, id="grizzly", slsp="missing")

# Regularise the time series
gbDateRef <-  as.POSIXct(format(grizzlyDate[1], "%Y-%m-%d 04:00"), tz="UTC")
cMov <- setNA(cMov, date.ref=gbDateRef, dt=4, units="hour")
cMov <- sett0(cMov, date.ref=gbDateRef, dt=4, units="hour")

#####
# Analyse the data
cRes <- movLikelihoods(cMov, PRdetails=TRUE, TAc=10, conts=FALSE, hspo=TRUE)
# CCRW_A, CCRW_L, TLW, BW, CRW
cAICc <- c(cRes$mleMov$CCRW["AICc"], cRes$mleMov$HSMMpo["AICc"], 
  cRes$mleMov$TLW["AICc"],
  cRes$mleMov$BW["AICc"], cRes$mleMov$CRW["AICc"])
# Delta AICc as in Table S2.2
cAICc - min(cAICc)

# p-value for test of absolute fit, based on Step length only (as in Table S2.2)
cRes$pseudoRes$Z["pval","SL_HSMMpo"]
