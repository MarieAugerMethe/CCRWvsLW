# Run the CCRW simulations of main manuscript

library(CCRWvsLW)

set.seed(123)

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

AICcRes <- matrix(ncol=7,nrow=nrow(parSeq)*nSim)
colnames(AICcRes) <- c("CCRW","LW","TLW","BW","CRW","TBW","TCRW")

prRes <- matrix(ncol=7,nrow=nrow(parSeq)*nSim)
colnames(prRes) <- c("CCRW","LW","TLW","BW","CRW","TBW","TCRW")

prResDet <- matrix(ncol=8,nrow=nrow(parSeq)*nSim)
colnames(prResDet) <- c("SL_CCRW", "SL_LW", "SL_TLW", "SL_E", "SL_TE", "TA_CCRW","TA_U", "TA_VM")

parRes <- matrix(ncol=7,nrow=nrow(parSeq)*nSim)
colnames(parRes) <- c("gII", "gEE", "lI", "lE", "kE", "dI", "a")

gIIciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(gIIciRes) <- c("simgII", "gII", "L95CI", "U95CI")

gEEciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(gEEciRes) <- c("simgEE", "gEE", "L95CI", "U95CI")

lIciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(lIciRes) <- c("simlI", "lI", "L95CI", "U95CI")

lEciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(lEciRes) <- c("simlE", "lE", "L95CI", "U95CI")

kEciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(kEciRes) <- c("simkE", "kE", "L95CI", "U95CI")

t0 <- proc.time()
for(i in 1:nrow(parSeq)){
#for(i in 1:5){
  for(j in 1:nSim){
  #for(j in 1:1){
    mov <- simmCCRW(500,parSeq[i,2], parSeq[i,3], lI, parSeq[i,1], parSeq[i,4], a, dI)
    movRes <- movLikelihoods(mov, PRdetails=TRUE, graph=FALSE)
    AICcRes[((i-1)*nSim+j),] <- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] 
    prRes[((i-1)*nSim+j),] <- movRes$pseudoRes$PR["pval",]
    prResDet[((i-1)*nSim+j),] <- movRes$pseudoRes$Z["pval",]
    parRes[((i-1)*nSim+j),] <- movRes$mle$CCRW[c(1:6,8)]
    gIIciRes[((i-1)*nSim+j),1] <- parSeq[i,2]
    gIIciRes[((i-1)*nSim+j),2:4] <- movRes$CI$CCRW[1,]
    gEEciRes[((i-1)*nSim+j),1] <- parSeq[i,3]
    gEEciRes[((i-1)*nSim+j),2:4] <- movRes$CI$CCRW[2,]
    lIciRes[((i-1)*nSim+j),1] <- lI
    lIciRes[((i-1)*nSim+j),2:4] <- movRes$CI$CCRW[3,]
    lEciRes[((i-1)*nSim+j),1] <- parSeq[i,1]
    lEciRes[((i-1)*nSim+j),2:4] <- movRes$CI$CCRW[4,]
    kEciRes[((i-1)*nSim+j),1] <- parSeq[i,4]
    kEciRes[((i-1)*nSim+j),2:4] <- movRes$CI$CCRW[5,]
  }
}
proc.time() - t0

write.csv(AICcRes, "CCRW_AICc_ResPlain.csv", row.names=FALSE)
write.csv(prRes, "CCRW_pr_ResPlain.csv", row.names=FALSE)
write.csv(prResDet, "CCRW_prDet_ResPlain.csv", row.names=FALSE)
write.csv(parRes, "CCRW_par_ResPlain.csv", row.names=FALSE)
write.csv(gIIciRes, "CCRW_gIIci_ResPlain.csv", row.names=FALSE)
write.csv(gEEciRes, "CCRW_gEEci_ResPlain.csv", row.names=FALSE)
write.csv(lIciRes, "CCRW_lIci_ResPlain.csv", row.names=FALSE)
write.csv(lEciRes, "CCRW_lEci_ResPlain.csv", row.names=FALSE)
write.csv(kEciRes, "CCRW_kEci_ResPlain.csv", row.names=FALSE)

library(plyr)
# How many NA due to underflow
warnings()
nrow(AICcRes[is.na(AICcRes[,1]),])/ nrow(AICcRes) #~2%
AICcMin <- apply(AICcRes, 1, which.min) # If it has an NA it just ignores that column for that simulation
AICcCount <- count(colnames(AICcRes)[AICcMin])
AICcCount

###
# Remove NA & focuss on CCRW, TLW, BW, CRW
naCCRW <- is.na(AICcRes[,"CCRW"])
AICcResMM <- AICcRes[,c("CCRW", "TLW", "BW", "CRW")]
AICcMinMM <- apply(AICcResMM, 1, min)
# Akaike weigth
expDAICcMM <- exp(-0.5*(AICcResMM - AICcMinMM))
Aw <- expDAICcMM/rowSums(expDAICcMM)

# Look at Akaike weigth (exclude simulations with NAs)
sum(Aw[!naCCRW ,"CCRW"] > 0.99)/nrow(Aw[!naCCRW,])*100
# Just checking if those lower results are due to step length being the same, at least partly!
#sum(Aw[rep(parSeq[,1],each=nSim) != lI,"CCRW"] > 0.99)/nrow(Aw[rep(parSeq[,1],each=nSim) != lI,])*100
sum(Aw[!naCCRW,"TLW"] > 0.01)/nrow(Aw[!naCCRW,])*100
sum(Aw[!naCCRW,"BW"] + Aw[!naCCRW,"CRW"] > 0.5)/nrow(Aw[!naCCRW,])*100

misExp <- cbind(lE=rep(parSeq[,1],each=nSim),
                gII=rep(parSeq[,2],each=nSim),
                gEE=rep(parSeq[,3],each=nSim),
                kE=rep(parSeq[,4],each=nSim),
         Aw[,"BW"] + Aw[,"CRW"] > 0.5)

count(misExp[misExp[,5]==1 & !naCCRW,1]) #lE - only for 0.005 and 0.01
count(misExp[misExp[,5]==1 & !naCCRW,2]) #gII - similarish across all par val
count(misExp[misExp[,5]==1 & !naCCRW,3]) #gEE - all par value affected but increase with lower values
count(misExp[misExp[,5]==1 & !naCCRW,4]) #kE - all par value affected but big increase with lower values

#########
# Look at whether the test of absolute fit work for TLW
sum(prRes[,"CCRW"] < 0.05)/nrow(prRes)
# This is huge, and thus we are rejecting the "true"/simulated model too often.
# This is likely partly driven by the turning angle distribution
sum(prResDet[,"SL_CCRW"] < 0.05)/nrow(prResDet)
# This is much lower, but still rejecting at a higher rate than 0.05
# This is partly because using the threshold angle change the path, 
# could suggests for this to use a more conservative threshold angle,
# if you don't want to reject the models (such as alpha 0.01)
sum(prResDet[,"SL_CCRW"] < 0.01)/nrow(prResDet)

#############
# Look at parameter estimates

layout(matrix(1:7,nrow=1))
par(las=1)
boxplot(parRes[,"gII"]~rep(parSeq[,"gII"],each=nSim))
boxplot(parRes[,"gEE"]~rep(parSeq[,"gEE"],each=nSim))
boxplot(parRes[,"lI"]~rep(lI,nrow(parRes)))
boxplot(parRes[,"lE"]~rep(parSeq[,"lE"],each=nSim))
boxplot(parRes[,"kE"]~rep(parSeq[,"kE"],each=nSim))
boxplot(parRes[,"a"]~rep(a,nrow(parRes)))
boxplot(parRes[,"dI"]~rep(dI,nrow(parRes)))

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
