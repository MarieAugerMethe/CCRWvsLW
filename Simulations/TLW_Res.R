# Analysis the results from TLW simulations
# You need to run TLW_loop before you can run this script

library(plyr)

##################################
# Read loop results files
AICcRes <- read.csv("TLW_AICc_ResPlain.csv")
prRes <- read.csv("TLW_pr_ResPlain.csv")
prResDet <- read.csv("TLW_prDet_ResPlain.csv")
parRes <- read.csv("TLW_par_ResPlain.csv")
ciRes <- read.csv("TLW_ci_ResPlain.csv")

################################
# Parameter used in the simulations
nSim <- 50
mu <- seq(1.1, 2.9, 0.1)
a <- 1
b <- c(100,1000,10000)
parSeq <- cbind(mu=rep(mu,length(b)), b=rep(b,each=length(mu)))

############################
# How many NA due to underflow
nrow(AICcRes[is.na(AICcRes[,1]),])/nrow(AICcRes)*100 # 107/2850 - so 3.75% (for CCRW TAc: 2/36000, so total: 0.28%) 
notNAccrw <- !is.na(AICcRes[,1])

############################
# Look at main models
AICcResMM <- AICcRes[,c("CCRW", "TLW", "BW", "CRW")]
AICcMinMM <- apply(AICcResMM, 1, min, na.rm=TRUE)
# Akaike weigth
expDAICcMM <- exp(-0.5*(AICcResMM - AICcMinMM))
Aw <- expDAICcMM/rowSums(expDAICcMM)

# Look at Akaike weigth (exclude simulations with NAs)
sum(Aw[notNAccrw,"TLW"] > 0.99)/nrow(Aw[notNAccrw,])*100
# Is there any case for misidentification
sum(Aw[notNAccrw,"CCRW"] > 0.01)/nrow(Aw[notNAccrw,])*100 # Small % of misidentification
sum(Aw[notNAccrw,"CCRW"] > 0.5)/nrow(Aw[notNAccrw,])*100 # Small % of misidentification, but not that bad

# Identify when the CCRW is misidentified for the CRW and BW
misExp <- cbind(mu=rep(parSeq[,1],each=nSim),
                b=rep(parSeq[,2],each=nSim),
                Aw[,"CCRW"] > 0.5)#0.1) 
misExp <- misExp[notNAccrw,]

count(misExp[misExp[,3]==1, 1]) #mu - 1.6, 2.6, 2.8
count(misExp[misExp[,3]==1, 2]) #b - mainly for short truncation points

##########
# Look at 7 models
AICcMin <- apply(AICcRes, 1, min, na.rm=TRUE)
# Akaike weigth
expDAICc <- exp(-0.5*(AICcRes - AICcMin))
Aw7 <- expDAICc/rowSums(expDAICc)

# Look at Akaike weigth (exclude simulations with NAs)

sum(Aw7[notNAccrw,"TLW"] + Aw7[notNAccrw,"LW"] > 0.99)/nrow(Aw7[notNAccrw,])*100 
sum(Aw7[notNAccrw,"CCRW"] > 0.01)/nrow(Aw7[notNAccrw,])*100
sum(Aw7[notNAccrw,"CCRW"] > 0.5)/nrow(Aw7[notNAccrw,])*100


#########
# Look at whether the test of absolute fit work for TLW
sum(prRes[,"TLW"] < 0.05)/nrow(prRes)

##############
# Look at quadractic CI
ciNaN <- is.nan(ciRes[,3]) | is.nan(ciRes[,4]) | is.na(ciRes[,3]) | is.na(ciRes[,4])
sum(ciRes[,3] <= ciRes[,1] & ciRes[,4] >= ciRes[,1], na.rm=TRUE)/nrow(ciRes)
sum(ciRes[,3] <= ciRes[,1] & ciRes[,4] >= ciRes[,1], na.rm=TRUE)/nrow(ciRes[!ciNaN,])
