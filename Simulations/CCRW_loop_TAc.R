# Run the CCRW simulations of Appendix with the local turn technique

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
    movRes <- movLikelihoods(mov, PRdetails=TRUE, graph=FALSE, TAc=10)
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

write.csv(AICcRes, "CCRW_AICc_ResS.csv", row.names=FALSE)
write.csv(prRes, "CCRW_pr_ResS.csv", row.names=FALSE)
write.csv(prResDet, "CCRW_prDet_ResS.csv", row.names=FALSE)
write.csv(parRes, "CCRW_par_ResS.csv", row.names=FALSE)
write.csv(gIIciRes, "CCRW_gIIci_ResS.csv", row.names=FALSE)
write.csv(gEEciRes, "CCRW_gEEci_ResS.csv", row.names=FALSE)
write.csv(lIciRes, "CCRW_lIci_ResS.csv", row.names=FALSE)
write.csv(lEciRes, "CCRW_lEci_ResS.csv", row.names=FALSE)
write.csv(kEciRes, "CCRW_kEci_ResS.csv", row.names=FALSE)

warnings()
