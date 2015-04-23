# Analysis the results from CCRW simulations of Appendix with the local turn technique
# You need to run CCRW_loop_TAC before you can run this script

library(plyr)

##################################
# Read loop results files
AICcRes <- read.csv("CCRW_AICc_Res.csv")
prRes <- read.csv("CCRW_pr_Res.csv")
prResDet <- read.csv("CCRW_prDet_Res.csv")
parRes <- read.csv("CCRW_par_Res.csv")
gIIciRes <- read.csv("CCRW_gIIci_Res.csv")
gEEciRes <- read.csv("CCRW_gEEci_Res.csv")
lIciRes <- read.csv("CCRW_lIci_Res.csv")
lEciRes <- read.csv("CCRW_lEci_Res.csv")
kEciRes <- read.csv("CCRW_kEci_Res.csv")

################################
# Parameter used in the simulations
nSim <- 50
dI <- 0.5
a <- 1
lI <- 0.01
lE <- c(0.01, 0.005,0.001,0.0005,0.0001)
gII <- c(0.6,0.7,0.8,0.9)
gEE <- seq(0.1,0.9,0.1)
kE <- c(0.5,1,5,10)
parSeq <- cbind(lE=rep(lE,each=length(gII)*length(gEE)*length(kE)), gII=rep(gII,each=length(gEE)*length(kE)),
                gEE=rep(gEE,each=length(kE)), kE=kE)

############################
# How many NA due to underflow
nrow(AICcRes[is.na(AICcRes[,1]),])/nrow(AICcRes) # 0/36000

############################
# Look at main models
AICcResMM <- AICcRes[,c("CCRW", "TLW", "BW", "CRW")]
AICcMinMM <- apply(AICcResMM, 1, min)
# Akaike weigth
expDAICcMM <- exp(-0.5*(AICcResMM - AICcMinMM))
Aw <- expDAICcMM/rowSums(expDAICcMM)

# Look at Akaike weigth
sum(Aw[,"CCRW"] > 0.99)/nrow(Aw)*100
# Is there any case for misidentification
sum(Aw[,"TLW"] > 0.01)/nrow(Aw)*100 # No misidentification with the TLW
sum(Aw[,"BW"] + Aw[,"CRW"] > 0.5)/nrow(Aw)*100 #Lots of misidentification with BW and CRW

# Identify when the CCRW is misidentified for the CRW and BW
misExp <- cbind(lE=rep(parSeq[,1],each=nSim),
                gII=rep(parSeq[,2],each=nSim),
                gEE=rep(parSeq[,3],each=nSim),
                kE=rep(parSeq[,4],each=nSim),
                Aw[,"BW"] + Aw[,"CRW"] > 0.5)

count(misExp[misExp[,5]==1, 1]) #lE
count(misExp[misExp[,5]==1, 2]) #gII
count(misExp[misExp[,5]==1, 3]) #gEE
count(misExp[misExp[,5]==1, 4]) #kE

#########
# Look at whether the test of absolute fit work for TLW
sum(prRes[,"CCRW"] < 0.05)/nrow(prRes)
# This is huge, and thus we are rejecting the "true"/simulated model too often.
# This is likely partly driven by the turning angle distribution
sum(prResDet[,"SL_CCRW"] < 0.05)/nrow(prResDet)
# This is much lower, but still rejecting at a higher rate than 0.05
# If you don't want to reject the models use a higher alpha (e.g., 0.01)
sum(prResDet[,"SL_CCRW"] < 0.01)/nrow(prResDet)

##############
# Look at quadractic CI

gIIciNaN <- is.nan(gIIciRes[,3]) | is.nan(gIIciRes[,4])
sum(gIIciRes[,3] <= gIIciRes[,1] & gIIciRes[,4] >= gIIciRes[,1], na.rm=TRUE)/nrow(gIIciRes)
sum(gIIciRes[,3] <= gIIciRes[,1] & gIIciRes[,4] >= gIIciRes[,1], na.rm=TRUE)/nrow(gIIciRes[!gIIciNaN,])

gEEciNaN <- is.nan(gEEciRes[,3]) | is.nan(gEEciRes[,4])
sum(gEEciRes[,3] <= gEEciRes[,1] & gEEciRes[,4] >= gEEciRes[,1], na.rm=TRUE)/nrow(gEEciRes)
sum(gEEciRes[,3] <= gEEciRes[,1] & gEEciRes[,4] >= gEEciRes[,1], na.rm=TRUE)/nrow(gEEciRes[!gEEciNaN,])

lIciNaN <- is.nan(lIciRes[,3]) | is.nan(lIciRes[,4])
sum(lIciRes[,3] <= lIciRes[,1] & lIciRes[,4] >= lIciRes[,1], na.rm=TRUE)/nrow(lIciRes)
sum(lIciRes[,3] <= lIciRes[,1] & lIciRes[,4] >= lIciRes[,1], na.rm=TRUE)/nrow(lIciRes[!lIciNaN,])

lEciNaN <- is.nan(lEciRes[,3]) | is.nan(lEciRes[,4])
sum(lEciRes[,3] <= lEciRes[,1] & lEciRes[,4] >= lEciRes[,1], na.rm=TRUE)/nrow(lEciRes)
sum(lEciRes[,3] <= lEciRes[,1] & lEciRes[,4] >= lEciRes[,1], na.rm=TRUE)/nrow(lEciRes[!lEciNaN,])

kEciNaN <- is.nan(kEciRes[,3]) | is.nan(kEciRes[,4])
sum(kEciRes[,3] <= kEciRes[,1] & kEciRes[,4] >= kEciRes[,1], na.rm=TRUE)/nrow(kEciRes)
sum(kEciRes[,3] <= kEciRes[,1] & kEciRes[,4] >= kEciRes[,1], na.rm=TRUE)/nrow(kEciRes[!kEciNaN,])

