fitGraph<- function(SL,TA,SLmin,SLmax,n,mleMov,glog="xy"){
	
	# PDF for step length of the models used for graphs and for G-test

	SLProbCCRW <- function(SL){
	# This is actually a CCRW_IM based in the stationary distribution calculated with CCRW
		mleMov$CCRW['I*'] * dexp((SL-SLmin), mleMov$CCRW[3]) +
			mleMov$CCRW['E*'] * dexp((SL-SLmin), mleMov$CCRW[4]) 
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


	TAProbCCRW <- function(TA_f){
		mleMov$CCRW['I*'] * dvm(TA_f, 0, 0) +
			mleMov$CCRW['E*'] * dvm(TA_f, 0, mleMov$CCRW[5]) 
	}

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
	# Need to add the if statement for no starting values giving good parameter values?
	if(!any(is.na(mleMov$CCRW[c('I*','E*','lI','lE')]))){
	  curve(SLProbCCRW, from=SLmin, to=SLmax, log=glog, add=TRUE, col=colc[5], lwd=2, lty=1, xpd=FALSE) 
	}
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
	if(!any(is.na(mleMov$CCRW[c('I*','E*','kE')]))){
	  curve(TAProbCCRW, from=-pi, to=pi, add=TRUE, col=colc[5], lwd=2, lty=1)
	}
	legend("topright",c("CCRW", "CRW/TCRW", "LW/TLW/","BW/TBW"),
	       col=c(colc[c(5,3,1)],"white"), lty=1, lwd=2, bty="n", seg.len=0.95)
  box()

	
  
}
