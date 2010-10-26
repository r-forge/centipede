plotNegBinPars <-
function(NegBinPar,mixP,x=(0:100))
  {
    B0 <- plogis(NegBinPar$logitB0)
    A0 <- exp(NegBinPar$logA0)
    B1 <- plogis(NegBinPar$logitB1)
    A1 <- exp(NegBinPar$logA1)
    Mean0 <- exp(NegBinPar$logA0)*exp(-NegBinPar$logitB0)
    Var0 <- Mean0/plogis(-NegBinPar$logitB0)
    Mean1 <- exp(NegBinPar$logA1)*exp(-NegBinPar$logitB1)
    Var1 <- Mean1/plogis(-NegBinPar$logitB1)
    lines(x,dnbinom(x,A0,B0),t='l',col='darkred',ylab="PMF",xlab="Num. Reads")
    lines(x,dnbinom(x,A1,B1),col='darkblue')
    if(!missing(mixP)){
      lines(x,(1-mixP)*dnbinom(x,A0,B0)+mixP*dnbinom(x,A1,B1),col='darkgreen')
    }
    legend("topright",c(paste("Not Bound:",round(Mean0,digits=2),round(Var0,digits=2)),paste("Bound:",round(Mean1,digits=2),round(Var1,digits=2))),col=c("darkred","darkblue"),lty=1)
  }

