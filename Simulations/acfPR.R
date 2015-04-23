# To look at the autocorelation in the pseudo-residuals as in the appendix

library(CCRWvsLW)

set.seed(782)
movltraj <- simmCCRW(500, 0.9, 0.9, 0.01,0.0001,10,1,0.5)
movRes <- movLikelihoods(movltraj, graph=FALSE, PRdetails=TRUE)

# To get at the actual pseudoresiduals
movD <- movFormat(movltraj)
PRres <- pseudo(movD$SL,movD$TA_C,movD$TA,movD$SLmin,movD$SLmax,
       movD$missL,movD$notMisLoc,movD$n,movRes$mle,PRdetails=FALSE,graph=TRUE,Uinfo=TRUE)
head(PRres$U)


# So to look at ACF for each of the main model
layout(matrix(1:8,nrow=4, byrow=TRUE))
par(las=1, mgp=c(2,0.5,0), tck=-0.03, mar=c(4,4,1.5,0.5))
sameAxis <- range(acf(PRres$U[,"U_SL_CCRW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_TA_CCRW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_TLW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_E"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_U"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_VM"],plot=FALSE)[[1]][-1])
acf(PRres$U[,"U_SL_CCRW"], main="", ylim=sameAxis, xlim=c(1,10))
title(main = "CCRW SL", line=0.3)
acf(PRres$U[,"U_TA_CCRW"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CCRW TA", line=0.3)
acf(PRres$U[,"U_TLW"], main = "", ylim=sameAxis, xlim=c(1,10))
title("TLW SL", line=0.3)
acf(PRres$U[,"U_U"], main = "", ylim=sameAxis, xlim=c(1,10))
title("TLW TA", line=0.3)
acf(PRres$U[,"U_E"], main="", ylim=sameAxis, xlim=c(1,10))
title(main = "BW SL", line=0.3)
acf(PRres$U[,"U_U"], main = "", ylim=sameAxis, xlim=c(1,10))
title("BW TA", line=0.3)
acf(PRres$U[,"U_E"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CRW SL", line=0.3)
acf(PRres$U[,"U_VM"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CRW TA", line=0.3)

# look at the impact of TAc

movRes <- movLikelihoods(movltraj, graph=FALSE, PRdetails=TRUE, TAc=10)

# To get at the actual pseudoresiduals
movD <- movFormat(movltraj, TAc=10)
PRres <- pseudo(movD$SL,movD$TA_C,movD$TA,movD$SLmin,movD$SLmax,
                movD$missL,movD$notMisLoc,movD$n,movRes$mle,PRdetails=FALSE,graph=TRUE,Uinfo=TRUE)
head(PRres$U)

# So to look at ACF for each of the main model
layout(matrix(1:8,nrow=4, byrow=TRUE))
par(las=1, mgp=c(2,0.5,0), tck=-0.03, mar=c(4,4,1.5,0.5))
sameAxis <- range(acf(PRres$U[,"U_SL_CCRW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_TA_CCRW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_TLW"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_E"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_U"],plot=FALSE)[[1]][-1],
                  acf(PRres$U[,"U_VM"],plot=FALSE)[[1]][-1])
acf(PRres$U[,"U_SL_CCRW"], main="", ylim=sameAxis, xlim=c(1,10))
title(main = "CCRW SL", line=0.3)
acf(PRres$U[,"U_TA_CCRW"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CCRW TA", line=0.3)
acf(PRres$U[,"U_TLW"], main = "", ylim=sameAxis, xlim=c(1,10))
title("TLW SL", line=0.3)
acf(PRres$U[,"U_U"], main = "", ylim=sameAxis, xlim=c(1,10))
title("TLW TA", line=0.3)
acf(PRres$U[,"U_E"], main="", ylim=sameAxis, xlim=c(1,10))
title(main = "BW SL", line=0.3)
acf(PRres$U[,"U_U"], main = "", ylim=sameAxis, xlim=c(1,10))
title("BW TA", line=0.3)
acf(PRres$U[,"U_E"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CRW SL", line=0.3)
acf(PRres$U[,"U_VM"], main = "", ylim=sameAxis, xlim=c(1,10))
title("CRW TA", line=0.3)
