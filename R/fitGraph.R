fitGraph<- function(SL, TA, SLmin, SLmax, n, mleMov, ccrwBm, glog="xy"){
	
  m <- c(10,10)
	# PDF for step length of the models used for graphs and for G-test
  if(ccrwBm == "o" | ccrwBm == "dn"){
    SLProbCCRW <- function(SL){
      # This is actually a CCRW_IM based in the stationary distribution calculated with CCRW
      mleMov$CCRW['I*'] * dexp((SL-SLmin), mleMov$CCRW[3]) +
        mleMov$CCRW['E*'] * dexp((SL-SLmin), mleMov$CCRW[4])
    }
    TAProbCCRW <- function(TA){
      mleMov$CCRW['I*'] * dvm(TA, 0, 0) +
        mleMov$CCRW['E*']* dvm(TA, 0, mleMov$CCRW["kE"])
    }
  }else if(ccrwBm == "ww"){
    SLProbCCRW <- function(SL){
      # This is actually a CCRW_IM based in the stationary distribution calculated with CCRW
      mleMov$CCRWww['I*'] * dweibull(SL,mleMov$CCRWww[5],mleMov$CCRWww[3]) +
        mleMov$CCRWww['E*'] * dweibull(SL,mleMov$CCRWww[6],mleMov$CCRWww[4])
    }
    TAProbCCRW <- function(TA){
      mleMov$CCRWww['I*'] * dwrcpauchy(TA, 0, 0) +
        mleMov$CCRWww['E*']* dwrpcauchy(TA, 0, mleMov$CCRWww["rE"])
    }
  }else if(ccrwBm == "hs"){
    gamma <- gen.Gamma.repar(m,mleMov$HSMM[c("gSI","gSE")],mleMov$HSMM[c("gPI","gPE")]) # Creating transition probility matrix
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))
    delta2 <- c(sum(delta[1:m[1]]),sum(delta[(m[1]+1):(m[2]+m[1])]))
    SLProbCCRW <- function(SL){
      delta2[1] * dweibull(SL,mleMov$HSMM["shI"],mleMov$HSMM["scI"]) +
        delta2[2] * dweibull(SL,mleMov$HSMM["shE"],mleMov$HSMM["scE"])
    }
    TAProbCCRW <- function(TA){
      delta2[1] * dwrcpauchy(TA, 0, 0) +
        delta2[2] * dwrpcauchy(TA, 0, mleMov$HSMM["rE"])
    }
  }else if(ccrwBm == "hsl"){
    gamma <- gen.Gamma.repar(m,mleMov$HSMMl[c("gSI","gSE")],mleMov$HSMMl[c("gP","gP")]) # Creating transition probility matrix
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))
    delta2 <- c(sum(delta[1:m[1]]),sum(delta[(m[1]+1):(m[2]+m[1])]))
    SLProbCCRW <- function(SL){
      delta2[1] * dweibull(SL,mleMov$HSMMl["shI"],mleMov$HSMMl["scI"]) +
               delta2[2] * dweibull(SL,mleMov$HSMMl["shE"],mleMov$HSMMl["scE"])
    }
    TAProbCCRW <- function(TA){
      delta2[1] * dwrcpauchy(TA, 0, 0) +
        delta2[2] * dwrpcauchy(TA, 0, mleMov$HSMMl["rE"])
    }
  }else if(ccrwBm == "hsp"){
    gamma <- gen.Gamma.pois(m,c(mleMov$HSMMp["laI"],mleMov$HSMMp["laE"])) # Creating transition probility matrix
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    delta2 <- c(sum(delta[1:m[1]]),sum(delta[(m[1]+1):(m[2]+m[1])]))
    SLProbCCRW <- function(SL){
      delta2[1] * dweibull(SL,mleMov$HSMMp["shI"],mleMov$HSMMp["scI"]) +
        delta2[2] * dweibull(SL,mleMov$HSMMp["shE"],mleMov$HSMMp["scE"])
    }
    TAProbCCRW <- function(TA){
      delta2[1] * dwrpcauchy(TA, 0, 0) +
        delta2[2] * dwrpcauchy(TA, 0, mleMov$HSMMp["rE"])
    }
  }else if(ccrwBm == "hspo"){
    gamma <- gen.Gamma.pois(m,c(mleMov$HSMMpo["laI"],mleMov$HSMMpo["laE"])) # Creating transition probility matrix
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    delta2 <- c(sum(delta[1:m[1]]),sum(delta[(m[1]+1):(m[2]+m[1])]))
    SLProbCCRW <- function(SL){
      delta2[1] * dexp((SL-SLmin), mleMov$HSMMpo["lI"]) +
               delta2[2] * dexp((SL-SLmin), mleMov$HSMMpo["lE"])
    }
    TAProbCCRW <- function(TA){
      delta2[1] * dvm(TA, 0, 0) +
        delta2[2] * dvm(TA, 0, mleMov$HSMMpo["kE"])
    }
  }
  

	
	SLProbLW <- function(SL){
		(mleMov$LW[1] - 1) * SLmin^(mleMov$LW[1]-1) * SL^(-mleMov$LW[1])
	}
		
	SLProbTLW <- function(SL){
		((mleMov$TLW[1] - 1)/(SLmin^(1-mleMov$TLW[1]) - SLmax^(1-mleMov$TLW[1])))* SL^(-mleMov$TLW[1])
	}
	
	# The BW and the CRW have the same PDF for the SL
	SLProbE <- function(SL){
		dexp((SL-SLmin),mleMov$BW[1])
	}

	# The TBW and the TCRW have the same PDF for the SL
	SLProbTE <- function(SL){
		(mleMov$TBW[1] * exp(-mleMov$TBW[1] *SL))/ (exp(-mleMov$TBW[1]*SLmin) - exp(-mleMov$TBW[1]*SLmax))
	}


	##
	# Turning angle

	# PDF for turning angle
	TAProbU <- function(TA_f){
		dvm(TA_f, 0, 0)
	}
	
	TAProbVM <- function(TA_f){
		dvm(TA_f, 0, mleMov$CRW[2])
	}

  ##
  # Colors
	colc <- c("grey77","grey77","grey48","grey48","black")
  
  #########################################################
  # Step lengh graph
  # Depends on whether you want log-log, or just log or not log at all
  
	par(mar=c(5.1,5.1,4.1,2.1), mgp=c(3.5,0.5,0),tcl=-0.3)
  if(glog=="" | glog=="y"){
    SLB <- seq(SLmin,SLmax,length.out=round(n/10))

  }else{
    SLB <- exp(seq(log(SLmin),log(SLmax),length.out=round(n/10)))
  }
	h <- hist(SL,breaks=SLB,plot=FALSE)
	suppressWarnings(plot(h$mids,h$density,log=glog,
	     xlab="Step length", ylab="PDF/Step Length density",
	     type="h", lwd=4, cex.lab=1.2, cex.axis=1.2, las=1, col="darkgrey", xpd=FALSE))
	curve(SLProbE, from=SLmin, to=SLmax,log=glog, add=TRUE, col=colc[1], lwd=2, lty=2, xpd=FALSE)
	curve(SLProbTE, from=SLmin, to=SLmax,log=glog, add=TRUE, col=colc[2], lwd=2, lty=1, xpd=FALSE)
	curve(SLProbLW, from=SLmin, to=SLmax,log=glog, add=TRUE, col=colc[3], lwd=2, lty=1, xpd=FALSE)
	curve(SLProbTLW, from=SLmin, to=SLmax,log=glog, add=TRUE, col=colc[4], lwd=2, lty=2, xpd=FALSE)
	tryCatch(curve(SLProbCCRW, from=SLmin, to=SLmax, log=glog, add=TRUE, col=colc[5], lwd=2, lty=1, xpd=FALSE),
	         error=function(e) warning("Cannot plot the best CCRW"))
	legend("topright",c("CCRW", "LW", "TLW", "BW/CRW", "TBW/TCRW"),
	       col=colc[5:1], lty=c(1,1,2,1,2), lwd=2, bty="n", seg.len=1.7)
  box()
  
	# TA
	h <- hist(TA,breaks=round(n/10),plot=FALSE)
	plot(h$mids,h$density,  xlim=c(-pi,pi), ylim=c(0,max(h$density)+max(h$density)/20),
	     xlab="Relative turning angle (radians)", ylab="PDF/Angle density",
	     type="h", lwd=4, cex.lab=1.2, cex.axis=1.2, las=1, col="darkgrey", xpd=FALSE)
	curve(TAProbU, from=-pi, to=pi, add=TRUE, col=colc[1], lwd=2, lty=1)
	curve(TAProbVM,  from=-pi, to=pi, add=TRUE, col=colc[3], lwd=2, lty=1)
	tryCatch(curve(TAProbCCRW, from=-pi, to=pi, add=TRUE, col=colc[5], lwd=2, lty=1),
	         error=function(e) warning("Cannot plot the best CCRW"))
	legend("topright",c("CCRW", "CRW/TCRW", "LW/TLW/","BW/TBW"),
	       col=c(colc[c(5,3,1)],"white"), lty=1, lwd=2, bty="n", seg.len=0.95)
  box()

	
  
}
