extractNegBinPar <-
function(NegBinPar)
  {
    Mean0 <- exp(NegBinPar$logA0)*exp(-NegBinPar$logitB0)
    Var0 <- Mean0/plogis(-NegBinPar$logitB0)
    Mean1 <- exp(NegBinPar$logA1)*exp(-NegBinPar$logitB1)
    Var1 <- Mean1/plogis(-NegBinPar$logitB1)
    list(Mean0=Mean0,Var0=Var0,Mean1=Mean1,Var1=Var1);
  }

