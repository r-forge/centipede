fLogit <-
function(BetaLogit,Y,Ez){
    ##myPi <- plogis(Y %*% BetaLogit)
    ## -sum(Ez*log(myPi)+(1-Ez)*log(1-myPi))
    ##-sum(Ez*(Y %*% BetaLogit)+ log(1-myPi))
    -sum(Ez*(Y %*% BetaLogit) - log(1+exp(Y %*% BetaLogit)))
  }

gLogit <-
function(BetaLogit,Y,Ez){
    myPi <- plogis(Y %*% BetaLogit)
    -t(Ez-myPi) %*% Y
  }

simpleLogitScore <-
function(X,Ez){
   BetaFit <- optim(c(0,0),fLogit,gLogit,Y=as.matrix(data.frame(IntCept=1,X=X)),Ez=Ez,method="BFGS",control=list(maxit=500),hessian=TRUE);
   logitSE <- sqrt(diag(solve(as.matrix(BetaFit$hessian))))
   Zlogit <- BetaFit$par[2]/logitSE[2]
   Zlogit
 }

