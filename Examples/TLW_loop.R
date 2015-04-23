# Run the TLW simulations of the main manuscript

library(CCRWvsLW)

set.seed(123)

nSim <- 50
mu <- seq(1.1, 2.9, 0.1)
a <- 1
b <- c(100,1000,10000)
parSeq <- cbind(rep(mu,length(b)), rep(b,each=length(mu)))

AICcRes <- matrix(ncol=7,nrow=nrow(parSeq)*nSim)
colnames(AICcRes) <- c("CCRW","LW","TLW","BW","CRW","TBW","TCRW")

prRes <- matrix(ncol=7,nrow=nrow(parSeq)*nSim)
colnames(prRes) <- c("CCRW","LW","TLW","BW","CRW","TBW","TCRW")

prResDet <- matrix(ncol=8,nrow=nrow(parSeq)*nSim)
colnames(prResDet) <- c("SL_CCRW", "SL_LW", "SL_TLW", "SL_E", "SL_TE", "TA_CCRW","TA_U", "TA_VM")

parRes <- matrix(ncol=3,nrow=nrow(parSeq)*nSim)
colnames(parRes) <- c("mu", "a", "b")

ciRes <- matrix(ncol=4,nrow=nrow(parSeq)*nSim)
colnames(ciRes) <- c("simMu", "mu", "L95CI", "U95CI")

t0 <- proc.time()
for(i in 1:nrow(parSeq)){
  for(j in 1:nSim){
    mov <- simmTLW(500,parSeq[i,1],a,parSeq[i,2])
    movRes <- movLikelihoods(mov, PRdetails=TRUE, graph=FALSE)
    AICcRes[((i-1)*nSim+j),] <- unlist(movRes$mle)[grep("AICc", names(unlist(movRes$mle)))] 
    prRes[((i-1)*nSim+j),] <- movRes$pseudoRes$PR["pval",]
    prResDet[((i-1)*nSim+j),] <- movRes$pseudoRes$Z["pval",]
    parRes[((i-1)*nSim+j),] <- movRes$mle$TLW[1:3]
    ciRes[((i-1)*nSim+j),1] <- parSeq[i,1]
    ciRes[((i-1)*nSim+j),2:4] <- movRes$CI$TLW
    print(paste("i=",i,",j=",j,sep=""))
  }
}
proc.time() - t0

write.csv(AICcRes, "TLW_AICc_ResPlain.csv", row.names=FALSE)
write.csv(prRes, "TLW_pr_ResPlain.csv", row.names=FALSE)
write.csv(prResDet, "TLW_prDet_ResPlain.csv", row.names=FALSE)
write.csv(parRes, "TLW_par_ResPlain.csv", row.names=FALSE)
write.csv(ciRes, "TLW_ci_ResPlain.csv", row.names=FALSE)

warnings()
