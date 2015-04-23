# Analysis the results from TLW simulations with local turn method
# You need to run TLW_loop_TAc before you can run this script

library(plyr)

##################################
# Read loop results files
AICcRes <- read.csv("TLW_AICc_Res.csv")
prRes <- read.csv("TLW_pr_Res.csv")
prResDet <- read.csv("TLW_prDet_Res.csv")
parRes <- read.csv("TLW_par_Res.csv")
ciRes <- read.csv("TLW_ci_Res.csv")

################################
# Parameter used in the simulations
nSim <- 50
mu <- seq(1.1, 2.9, 0.1)
a <- 1
b <- c(100,1000,10000)
parSeq <- cbind(mu=rep(mu,length(b)), b=rep(b,each=length(mu)))

############################
# How many NA due to underflow
nrow(AICcRes[is.na(AICcRes[,1]),])/nrow(AICcRes)*100 # 90/2850 - so 3.16% (for CCRW TAc: 0/36000, so total: 90/(2850+36000)) 
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
sum(Aw[notNAccrw,"CCRW"] > 0.01)/nrow(Aw[notNAccrw,])*100
sum(Aw[notNAccrw,"CCRW"] > 0.5)/nrow(Aw[notNAccrw,])*100

# Identify when the CCRW is misidentified for the CRW and BW
misExp <- cbind(mu=rep(parSeq[,1],each=nSim),
                b=rep(parSeq[,2],each=nSim),
                Aw[,"CCRW"] > 0.5)
misExp <- misExp[notNAccrw,]

count(misExp[misExp[,3]==1, 1]) #mu
count(misExp[misExp[,3]==1, 2]) #b

#########
# Look at whether the test of absolute fit work for TLW
sum(prRes[,"TLW"] < 0.05)/nrow(prRes)
# This is huge, and thus we are rejecting the "true"/simulated model too often.
# This is likely partly driven by the turning angle distribution
sum(prResDet[,"SL_TLW"] < 0.05)/nrow(prResDet)
# This is much lower, but still rejecting at a higher rate than 0.05
# This is partly because using the threshold angle change the path. 
sum(prResDet[,"SL_TLW"] < 0.01)/nrow(prResDet)
sum(prResDet[,"TA_U"] < 0.05)/nrow(prResDet)

##############
# Look at quadractic CI
ciNaN <- is.nan(ciRes[,3]) | is.nan(ciRes[,4]) | is.na(ciRes[,3]) | is.na(ciRes[,4])
sum(ciRes[,3] <= ciRes[,1] & ciRes[,4] >= ciRes[,1], na.rm=TRUE)/nrow(ciRes)
sum(ciRes[,3] <= ciRes[,1] & ciRes[,4] >= ciRes[,1], na.rm=TRUE)/nrow(ciRes[!ciNaN,])
