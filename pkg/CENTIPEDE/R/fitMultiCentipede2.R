## X data -- rows location , -- columns base site
## Y covariates -- rows location , -- column covariate value, add a one for intercept point.
## Lambda -- profile parameters initial values
## BetaLogit -- Covariates parameters initial values
##' @param Xlist List of data matrices
##' @param Y Covariate data matrix
##' @param LambdaParList 
##' @param BetaLogit
##' @param NegBinParList
##' @param sweeps
##' @param DampLambda
##' @param NegBinPiShrink

## TODO: LambdaParList --> Lambda!!

##' @return ...
fitMultiCentipede <- function(Xlist,Y, KmerSeq,Mappability, KmerCutPr, BetaLogit, LambdaParList, NegBinParList,sweeps=50,DampLambda=0.01,DampNegBin=0.0,TrimP=0.0001,NRiter=5) {
  LogLikEvo <- vector(length=sweeps)

  ## Checking Number of Matrices... and dimensionality match
  stopifnot(is.list(Xlist)); 
  NumTissues <- length(Xlist); 
  stopifnot(NumTissues>0);
   
  ## Y<-as.matrix(Y)
  stopifnot(sapply(Xlist,is.matrix)); #Either this or convert to matrices
    
  L <- dim(Xlist[[1]])[1]  #Number of rows has to be the same on all matrices
  stopifnot(sapply(Xlist,function(i) {dim(i)[1]==L}))

  ## TOREMOVE!!
  Slist <- sapply(Xlist,function(i) {dim(i)[2]})  
  #Number of columns cannot differ for the multi-tissue model
  S <- dim(Xlist[[1]])[2]  #Number of cols has to be the same on all matrices
  stopifnot(sapply(Xlist,function(i) {dim(i)[2]==S}))

  if(missing(Y)){
    Y <- matrix(1,L); ## Intercept point. 
  }
  stopifnot(is.matrix(Y));
  stopifnot(dim(Y)[1]==L)
  
  K<-dim(Y)[2]
     
  Rlist <- lapply(Xlist,rowSums) #Total number of reads
  Rlist <- lapply(Rlist,as.matrix)

  ##Sequence preference...
  if(missing(KmerSeq) && missing(KmerCutPr)){
    RhoMat=1;
    Pcut=1;
    SeqDepFlag=FALSE;
  }else{
    SeqDepFlag=TRUE;
    stopifnot(!missing(KmerSeq) && !missing(KmerCutPr));
    KmerSeq <- as.matrix(KmerSeq)
    RhoMat <- t(apply(KmerSeq,1,function(row){ll <- KmerCutPr[(row)+1]; ll}))
    Pcut <- rowSums(RhoMat)
    Pcut=1;
    RhoMat <- t(apply(KmerSeq,1,function(row){ll <- KmerCutPr[(row)+1]; ll/sum(ll)}))
    RhoMat <- RhoMat*(1-DampLambda)+DampLambda*1/S
    ##RhoMat <- RhoMat*(1-0.01)+0.01*1/S
  }

  if(missing(Mappability)){
    MappabilityFlag=FALSE;
  }else{
    MappabilityFlag=TRUE;
    RhoMat <- RhoMat*Mappability;
    ##RhoMat <- t(apply(RhoMat,1,function(row){row/sum(row)}))
    ##RhoMat <- RhoMat*(1-DampLambda)+DampLambda*1/S
    ##RhoMat <- RhoMat*(1-0.01)+0.01*1/S
  }

  
  
  ## Preparing internal functions for the optimizers.
  ## Preparing functions for the logistic fit using R's BGFS -- 
  ## Objective function
  fLogit <- function(BetaLogit,Y,Ez){
    ##myPi <- plogis(Y %*% BetaLogit)
    ## -sum(Ez*log(myPi)+(1-Ez)*log(1-myPi))
    ##-sum(Ez*(Y %*% BetaLogit)+ log(1-myPi))
    -sum(Ez*(Y %*% BetaLogit) - log(1+exp(Y %*% BetaLogit)))
  }
  ## Gradient 
  gLogit <- function(BetaLogit,Y,Ez){
    myPi <- plogis(Y %*% BetaLogit)
    -t(Ez-myPi) %*% Y
  }
  ## Preparing functions for the Negative Binomial alpha_0 using R's BGFS -- 
  ## Objective function
  ## fNegBin <- function(AandB,R,Ez){
  ##   A<-AandB[1]
  ##   B<-AandB[2]
  ##   -(sum( Ez * (lgamma(R+A) - lgamma(A) + A*log(B) + R*log(1-B)) ))
  ## }
  ## gNegBin <- function(AandB,R,Ez){
  ##   A<-AandB[1]
  ##   B<-AandB[2]
  ##   gA <- -(sum( Ez * (digamma(R+A) - digamma(A) + log(B)) ))
  ##   gB <- -(sum( Ez * (A/B - R/(1-B)) ))
  ##   c(gA,gB)
  ## }
  ## Alternative parametrization !!!
  ## Objective function
  fNegBin2 <- function(NegBinPar, R, Ez=1,Pcut=1) {
    A <- exp(NegBinPar[1])
    logitB <- NegBinPar[2]-log(Pcut);
    if ((A < 1e+200) && (A > 1e-200)) {
      LogLik <- (sum(Ez * (lgamma(R + A) - lgamma(A) + 
                            A * logitB - (R + A) * log(1 + exp(logitB)))))
    }
    else {
      LogLik <- -1e+300
    }
    -LogLik
  }
  ## Gradient function
  gNegBin2 <- function(NegBinPar, R, Ez=1,Pcut=1) {
    A <- exp(NegBinPar[1])
    logitB <- NegBinPar[2]-log(Pcut)
    B <- plogis(logitB)
    gLogA <- -(A * sum(Ez * (digamma(R + A) - digamma(A) + 
                             logitB - log(1 + exp(logitB)))))
    gLogitB <- -(sum(Ez * (A - (A + R) * B)))
    c(gLogA, gLogitB)
  }

  ## Preparing functions for the Multinomial with seq dependece BGFS -- 
  ## Objective function
  fMultSeqDep <- function(LogPar,Xs,EzRl,rho){
    lambda <- exp(LogPar);
    lambda <- lambda/sum(lambda);
    ##aux <- apply(rho,1,function(rhorow){sum(rhorow*lambda)});
    aux <- rho %*% lambda
    -(sum(log(lambda)*Xs)-sum(EzRl*log(aux)))
  }
  gMultSeqDep <- function(LogPar,Xs,EzRl,rho){
    lambda <- exp(LogPar);
    lambda <- lambda/sum(lambda);
    ##aux <- apply(rho,1,function(rhorow){sum(rhorow*lambda)});
    aux <- rho %*% lambda
    ##-(Xs-lambda*colSums(EzRl/aux*rho))
    -(Xs-lambda*t(t(EzRl/aux) %*% rho))
  }

  ################################################################
  ################################################################

  compLambda <- function(X,Ez,DampLambda=0.0){
    S <- dim(X)[2];
    if(S>1){
      ##    Ez <- Ez*(1-DampLambda)+0.5*(DampLambda); ## Dampening experiment
      Lambda <- t(t(Ez) %*% X);
      Lambda <- Lambda/sum(Lambda);
      Lambda <- Lambda*(1-DampLambda) + (DampLambda)*1/S;
    }else{
      Lambda <- 1
    }
    Lambda
  }
  
  compLambdaSeqDep <- function(LambdaIni,X,R,Ez,DampLambda=0.0,rho){
    S <- dim(X)[2];
    if(S>1){
      ##    Ez <- Ez*(1-DampLambda)+0.5*(DampLambda); ## Dampening experiment
      ##Xs <- t(t(Ez) %*% X);
      ##EzRl <- Ez*Rl;
      ##warning(LambdaIni);
      LogParIni <- log(LambdaIni)
      optPar <- optim(LogParIni,fMultSeqDep,gMultSeqDep,
                      Xs = t(t(Ez) %*% X),EzRl = Ez*R,rho = rho,
                      method="BFGS",control = list(maxit = NRiter))
      LogPar <- optPar$par
      Lambda <- exp(LogPar);
      ##sum(Lambda)
      Lambda <- Lambda/sum(Lambda);
      Lambda <- Lambda*(1-DampLambda) + (DampLambda)*1/S;
    }else{
      Lambda <- 1
    }
    as.matrix(Lambda)
  }

  ## Do I need this anymore????
  compLambdaNoMix <- function(X,DampLambda=0.0){
    S <- dim(X)[2];
    if(S>1){
      ##    Ez <- Ez*LambPi+0.5*(1-LambPi); ## Dampening experiment
      Lambda <- colSums(X);
      Lambda <- Lambda/sum(Lambda);
      Lambda <- Lambda*(1-DampLambda) + (DampLambda)*1/S;
    }else{
      Lambda <- 1
    }
    Lambda
  }

  ## Initialization of the Negative Binomial, looks a bit horrible...
  ## But it is simply method of moments plus some clipping to aboid extreme values. 
  initNegBinParamsMixture <- function(R,Ez,DampNegBin=0.0){
    Ez <- Ez*(1-DampNegBin)+0.5*(DampNegBin);
    myM<-sum(R*Ez)/sum(Ez)
    myV<-sum(R^2 * Ez)/sum(Ez) - myM^2
    B1<-myM/myV;
        if(B1 > 0.999999){
        B1 <- 0.999999   
        }      
    A1<-myM*B1/(1-B1);
    ##
    myM<-sum(R*(1-Ez))/sum(1-Ez)
    myV<-sum(R^2 * (1-Ez))/sum(1-Ez) - myM^2
    B0<-myM/myV;
	if(B0 > 0.999999){
	B0 <- 0.999999
	}
    A0<-myM*B0/(1-B0);
    ##
    list(logitB0=qlogis(B0),logA0=log(A0),logitB1=qlogis(B1),logA1=log(A1))
  }
  
  compNegBinParamsMixture <- function(R,NegBinPar,Ez,Pcut=1,DampNegBin=0.0,NRiter=5){
    ## Negative Binomial dampening experiment!
    if(DampNegBin<1){
      Ez <- Ez*(1-DampNegBin)+0.5*(DampNegBin);
      ## Fast update first
      NegBinPar$logitB1 <- NegBinPar$logA1 + log(sum(Ez)) - log(sum(R*Ez))
      NegBinPar$logitB0 <- NegBinPar$logA0 + log(sum(1-Ez)) - log(sum(R*(1-Ez)))
    
      NegBinPar0 <- c(NegBinPar$logA0,NegBinPar$logitB0)
      NegBinPar1 <- c(NegBinPar$logA1,NegBinPar$logitB1)
      ## M-Step NegBinom parameters 
      NegBinPar0 <- optim(NegBinPar0, fNegBin2, gNegBin2, R=R, Ez=(1-Ez), Pcut=Pcut, method="L-BFGS-B",lower=c(-20,-20),upper=c(20,20),control=list(maxit=NRiter))
      NegBinPar1 <- optim(NegBinPar1, fNegBin2, gNegBin2, R=R, Ez=Ez    , Pcut=Pcut, method="L-BFGS-B",lower=c(-20,-20),upper=c(20,20),control=list(maxit=NRiter))
    
      NegBinPar$logA0 <- NegBinPar0$par[1];
      NegBinPar$logitB0 <- NegBinPar0$par[2];
      NegBinPar$logA1 <- NegBinPar1$par[1];
      NegBinPar$logitB1 <- NegBinPar1$par[2];
    }
    list(NegBinPar)
  }

  
  ## Without mixutre, NRiter should be higher since I will optimize only once...
  compNegBinParams <- function(R,NegBinPar,NRiter=5){
    ##str(NegBinPar)
    if(missing(NegBinPar)){
      myM<-mean(R);
      myV<-var(R);
      B<-myM/myV;
      A<-myM*B/(1-B);
      NegBinPar <- c(log(A),qlogis(B));
    }else{
      NegBinPar <- c(NegBinPar$logA,NegBinPar$logitB)
    }
    NegBinPar <- optim(NegBinPar, fNegBin2, gNegBin2, R=R, Ez=1 , method="L-BFGS-B",lower=c(-20,-20),upper=c(20,20),control=list(maxit=NRiter*100))$par   
    ##str(NegBinPar)
    list(logitB=NegBinPar[2],logA=NegBinPar[1])
  }

  shrink <- function(x,T){
    x[abs(x)<=T] <- 0
    x[x<(-T)] <-  x[x<(-T)]+T
    x[x>( T)] <-  x[x>( T)]-T
    x
  }
  clipExtremes <- function(x,TrimP){
    Q <- quantile(x,c(TrimP,1-TrimP))
    x[x<Q[1]] <- Q[1]
    x[x>Q[2]] <- Q[2]
    x
  }

  ## Added sequence stuff here... THINK IF THERE IS A BETTER WAY
  compNegBinLogRatio <- function(R,NegBinPar,TrimP=0.0001,Pcut=1,DampNegBin=0.0){
    if(DampNegBin<1){
      A0 <- exp(NegBinPar$logA0)
      A1 <- exp(NegBinPar$logA1)
      logitB0 <- NegBinPar$logitB0-log(Pcut)
      logitB1 <- NegBinPar$logitB1-log(Pcut)
      NegBinLogRatio <-  lgamma(A0) - lgamma(R+A0) - lgamma(A1) + lgamma(R+A1) +
        A1*logitB1 - (R+A1)*log(1+exp(logitB1)) - A0*logitB0 + (R+A0)*log(1+exp(logitB0));
      NegBinLogRatio <- clipExtremes(NegBinLogRatio,TrimP)
    }else{
      NegBinLogRatio <- R*0.0;
    }
    NegBinLogRatio
  }

  compMultLogRatioOcc2 <- function(X,Lambda,occ=1){
    Lambda <- Lambda * occ + (1-occ) * (1/length(Lambda));
    MultiNomLogRatio <-  X %*% log(Lambda*length(Lambda)); ## Multinomial for 1 and 0
  }
  compMultLogRatioOccPmax2 <- function(X,Lambda){
    tt <- compMultLogRatioOcc2(X,Lambda,1.0)
    ##tt <- pmax(tt,compMultLogRatioOcc(X,Lambda,0.3))
    ##tt <- pmax(tt,compMultLogRatioOcc(X,Lambda,0.4))
    tt <- pmax(tt,compMultLogRatioOcc2(X,Lambda,0.5))
    tt <- pmax(tt,compMultLogRatioOcc2(X,Lambda,0.6))
    tt <- pmax(tt,compMultLogRatioOcc2(X,Lambda,0.7))
    tt <- pmax(tt,compMultLogRatioOcc2(X,Lambda,0.8))
    tt <- pmax(tt,compMultLogRatioOcc2(X,Lambda,0.9))
    tt
  }  
  compMultiNomLogRatio <- function(X,Lambda,TrimP=0.0001){
    S <- dim(X)[2];
    if(S>1){
      ##MultiNomLogRatio <-  X %*% log(Lambda*S); ## Multinomial for 1 and 0
      MultiNomLogRatio <- compMultLogRatioOccPmax2(X,Lambda)      
      MultiNomLogRatio <- clipExtremes(MultiNomLogRatio,TrimP)
    }else{
      MultiNomLogRatio <- matrix(0,dim(X)[1],1)
    }
    MultiNomLogRatio
  }

  compMultLogRatioOcc <- function(X,R,rho,Lambda,occ=1){
    Lambda <- Lambda * occ + (1-occ) * (1/length(Lambda));
    logaux <- log(rho %*% Lambda)
    MultiNomLogRatio <-  X %*% log(Lambda) - R*logaux; ## Multinomial for 1 and 0
  }
  compMultLogRatioOccPmax <- function(X,R,rho,Lambda){
    tt <- compMultLogRatioOcc(X,R,rho,Lambda,1.0)
    ##tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.3))
    ##tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.4))
    tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.5))
    tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.6))
    tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.7))
    tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.8))
    tt <- pmax(tt,compMultLogRatioOcc(X,R,rho,Lambda,0.9))
    tt
  }
  compMultiNomLogRatioSeqDep <- function(X,R,Lambda,rho,TrimP=0.0001){
    S <- dim(X)[2];
    if(S>1){
      logaux <- log(rho %*% Lambda)
      MultiNomLogRatio <-  X %*% log(Lambda) - R*logaux; ## Multinomial for 1 and 0
      ##
      ##MultiNomLogRatio <-  X %*% log(Lambda) - rowSums(X*log(rho));
      ##
      ##MultiNomLogRatio <- compMultLogRatioOccPmax(X,R,rho,Lambda)
      ##
      MultiNomLogRatio <- clipExtremes(MultiNomLogRatio,TrimP)
    }else{
      MultiNomLogRatio <- matrix(0,dim(X)[1],1)
    }
    MultiNomLogRatio
  }

  compLogLik0 <- function(R,NegBinPar,S,TrimP=0.0001,DampNegBin=0.0){
    if(DampNegBin<1){
      A0 <- exp(NegBinPar$logA0)
      logitB0 <- NegBinPar$logitB0
      LogLik0 <-  lgamma(R+A0) - lgamma(A0) + A0*logitB0 - (R+A0)*log(1+exp(logitB0)) - R * log(S);
      LogLik0 <- clipExtremes(LogLik0,TrimP)
    }else{
      LogLik0 <- - R * log(S);
    }
  }
  compLogLik0SeqDep <- function(X,R,NegBinPar,rho,TrimP=0.0001,DampNegBin=0.0){
    if(DampNegBin<1){
      A0 <- exp(NegBinPar$logA0)
      logitB0 <- NegBinPar$logitB0
      LogLik0 <-  lgamma(R+A0) - lgamma(A0) + A0*logitB0 - (R+A0)*log(1+exp(logitB0)) + rowSums(X * log(rho));
      LogLik0 <- clipExtremes(LogLik0,TrimP)
    }else{
      LogLik0 <- rowSums(X * log(rho));
    }
  }

  compLogLik1 <- function(X,Lambda,NegBinPar,TrimP=0.0001,DampNegBin=0.0){
    S <- dim(X)[2];
    R <- rowSums(X);
    if(S>1){
      MultiNomLogLik <-  X %*% log(Lambda); ## Multinomial for 1 and 0
      ## MultiNomLogLik <- clipExtremes(MultiNomLogRatio,TrimP)
    }else{
      MultiNomLogLik <- matrix(0,dim(X)[1],1)
    }
    if(DampNegBin<1){
      A <- exp(NegBinPar$logA)
      logitB <- NegBinPar$logitB
      NegBinLogLik <-  lgamma(R+A) - lgamma(A) + A*logitB - (R+A)*log(1+exp(logitB));
    }else{
      NegBinLogLik <- 0.0;
    }
    LogLik1 <- NegBinLogLik+MultiNomLogLik;
    
    clipExtremes(LogLik1,TrimP)
  }
  compLogLik1SeqDep <- function(X,Lambda,NegBinPar,rho,TrimP=0.0001,DampNegBin=0.0){
    S <- dim(X)[2];
    R <- rowSums(X);
    if(S>1){
      ##MultiNomLogLik <-  X %*% log(Lambda); ## Multinomial for 1 and 0
      MultiNomLogLik <- compMultiNomLogRatioSeqDep(X,R,Lambda,rho,TrimP) + rowSums(X * log(rho)) 
      ## MultiNomLogLik <- clipExtremes(MultiNomLogRatio,TrimP)
    }else{
      MultiNomLogLik <- matrix(0,dim(X)[1],1)
    }
    if(DampNegBin<1){
      A <- exp(NegBinPar$logA)
      logitB <- NegBinPar$logitB
      NegBinLogLik <-  lgamma(R+A) - lgamma(A) + A*logitB - (R+A)*log(1+exp(logitB));
    }else{
      NegBinLogLik <- 0.0;
    }
    LogLik1 <- NegBinLogLik+MultiNomLogLik;
    
    clipExtremes(LogLik1,TrimP)
  }  


################################################################
################################################################
## Initializiation NEW
################################################################  
  cat("Initialization of the parameters:\n")

  ## Could add all-tissues together...
  Ez <- (Rlist[[1]]> quantile(Rlist[[1]],0.9))+0.0 ## Assuming Dnase is the first one...

  ## BetaLogit initialization, ....  
  if(missing(BetaLogit)){    
    BetaLogit <- matrix(0,K,1)
    BetaFit <- optim(BetaLogit,fLogit,gLogit,Y=Y,Ez=Ez,method="BFGS");
    BetaLogit<-BetaFit$par;    
  }

  row.names(BetaLogit) <- colnames(Y);
  
  #PriorLogRatio <- Y %*% BetaLogit;
  #Ez<-plogis(PriorLogRatio);

  ## Lambda initialization...
  ## M-Step Lambda 
  if(missing(LambdaParList)){
    LambdaParList <- lapply(Xlist,compLambda,Ez=Ez,DampLambda=DampLambda);
  }
  
## NegBin initialization...   
  ## M-Step NegBinParameters

  if(missing(NegBinParList)){
    NegBinParList <- lapply(Rlist,initNegBinParamsMixture,Ez=Ez,DampNegBin=DampNegBin);
    NegBinParList <- mapply(compNegBinParamsMixture,Rlist,NegBinParList,MoreArgs=list(Ez=Ez,DampNegBin=DampNegBin));
  }

  ########################
  ## No mixture Model ##
  
  NegBinParNoMixList <- lapply(Rlist,compNegBinParams,NRiter=100)

  LambdaParNoMixList <- lapply(Xlist,compLambdaNoMix,DampLambda=DampLambda)

  LogLikNoMix <- sum(mapply(compLogLik1,Xlist,LambdaParNoMixList,NegBinParNoMixList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)))

  LogLikNoMixNoNegBin <- sum(mapply(compLogLik1,Rlist,LambdaParNoMixList,NegBinParNoMixList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)))
  
############################################################
  
 ## Since I will do the difference with the backgraound only model I don't need this anymore!...
 MyLogLikConst <- 0 #( -log(S)*sum(X)); # - sum(lgamma(1+X)) - log(S)*sum(X)  #These guys are constants, so if we take them out, we go faster...

 ## First E-step outside the loop
 PriorLogRatio <- Y %*% BetaLogit;
 if(SeqDepFlag){
   NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,Pcut=Pcut,DampNegBin=DampNegBin)));
   MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatioSeqDep,Xlist,Rlist,LambdaParList,MoreArgs=list(rho=RhoMat, TrimP=TrimP)));
 }else{
   NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)));
   MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatio,Xlist,LambdaParList,MoreArgs=list(TrimP=TrimP)));
 }
 LogRatios <- PriorLogRatio+NegBinLogRatio+MultiNomLogRatio;
 Ez<-plogis(LogRatios);

 LogLikEvo <- 0;
 LogLik0 <- sum(mapply(compLogLik0,Rlist,NegBinParList,Slist,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(Y %*% BetaLogit)));
  Aux <- log(1+exp(LogRatios));
  Aux[LogRatios>18] <- LogRatios[LogRatios>18]; ## makes the thing more stable if LogRatio very large! so exp(LogRatio) does not explode
  LogLikEvo <- sum(Aux) + LogLik0; ## + sum(log(1-myPi)) - L*lgamma(A0) + L*log(B0)*A0 + sum(lgamma(R+A0)) + log(1-B0)*sum(R) + MyLogLikConst 
##  - sum(R*log(S))
  
#########################################################################################################
  ##  EM loop starts here ##
#########################################################################################################
  cat("\n Starting EM loop:\n")
##  for(i in 1:sweeps) {
  IterateEM <- TRUE
  ConvergedEM <- FALSE
  i <- 0;
  if(sweeps==0){
    IterateEM <- FALSE;
    i <- 1;
  }
  
  while(IterateEM){
    cat("----------------------------------------------------------------\n");
    i <- i+1;
    cat(" EM iteration ",i,":\n")
    ## E-Step 
    
    
    ## M-Step Lambda
    OldLambdas <- LambdaParList;
    if(SeqDepFlag){
      ##        compLambdaSeqDep <- function(LambdaIni,X,R,Ez,DampLambda=0.0,rho){
      LambdaParList <- mapply(compLambdaSeqDep, LambdaParList, Xlist,Rlist,MoreArgs=list(Ez=Ez,DampLambda=DampLambda,rho=RhoMat),SIMPLIFY=FALSE);
      ##LambdaParList <- lapply(Xlist,compLambda,Ez=Ez,DampLambda=DampLambda);
    }else{
      LambdaParList <- lapply(Xlist,compLambda,Ez=Ez,DampLambda=DampLambda);
    }
    
    ## M-Step NegBin
    OldNegBinPars <- NegBinParList;
    ##NegBinParList <- lapply(Rlist,initNegBinParamsMixture,Ez=Ez);
    NegBinParList <- mapply(compNegBinParamsMixture,Rlist,NegBinParList,MoreArgs=list(Ez=Ez,NRiter=NRiter,DampNegBin=DampNegBin,Pcut=Pcut));

    ## M-step Logistic,
    OldBetaLogit <- BetaLogit;
##    BetaFit <- optim(BetaLogit,fLogit,gLogit,Y=Y,Ez=Ez,method="L-BFGS-B",lower=rep(-20,K),upper=rep(20,K),control=list(maxit=5));
    BetaFit <- optim(BetaLogit,fLogit,gLogit,Y=Y,Ez=Ez,method="BFGS",control=list(maxit=NRiter));
    BetaLogit<-BetaFit$par;

    ## Monitor InComplete LogLik changes
    myPi <- plogis(Y %*% BetaLogit)

    PriorLogRatio <- Y %*% BetaLogit;
    if(SeqDepFlag){
      NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,Pcut=Pcut,DampNegBin=DampNegBin)));
      ##browser();
      MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatioSeqDep,Xlist,Rlist,LambdaParList,MoreArgs=list(rho=RhoMat, TrimP=TrimP)));
    }else{
      NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)));
      MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatio,Xlist,LambdaParList,MoreArgs=list(TrimP=TrimP)));
    }
    LogRatios <- PriorLogRatio+NegBinLogRatio+MultiNomLogRatio;
    Ez<-plogis(LogRatios);

    if(SeqDepFlag){
      LogLik0 <- sum(mapply(compLogLik0,Rlist,NegBinParList,Slist,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(Y %*% BetaLogit)));
    }else{
      LogLik0 <- sum(mapply(compLogLik0SeqDep,Xlist,Rlist,NegBinParList,MoreArgs=list(rho=RhoMat,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(Y %*% BetaLogit)));
    }
    
    Aux <- log(1+exp(LogRatios));
    Aux[LogRatios>18] <- LogRatios[LogRatios>18]; ## makes the thing more stable if LogRatio very large! so exp(LogRatio) does not explode
    LogLikEvo[i] <- sum(Aux) + LogLik0; ## + sum(log(1-myPi)) - L*lgamma(A0) + L*log(B0)*A0 + sum(lgamma(R+A0)) + log(1-B0)*sum(R) + MyLogLikConst 
##  - sum(R*log(S))

    ## Parameter Convergence monitoring..... 
    LambdaChange <-sapply(mapply(function(i,j){abs(j-i)},OldLambdas,LambdaParList),max);
    NegBinChange <-as.matrix(mapply(function(i,j){(as.data.frame(j)-as.data.frame(i))},OldNegBinPars,NegBinParList))
    LogitChange <- (BetaLogit-OldBetaLogit)
    MaxParamChange <- c(max(abs(LambdaChange)),max(abs(unlist(NegBinChange))),max(abs(LogitChange)))
    str(MaxParamChange)

    print("Max Lambda Changes:");
    if(NumTissues==1){
      print(max(abs(LambdaChange)));
    }else{
      print(LambdaChange);
    }
    print("Neg Bin Changes:");
    print(NegBinChange);
    print("Neg Bin New Values:");
    print(as.data.frame(lapply(NegBinParList,cbind)))
    print("Logistic Changes:");
    print(LogitChange);
    print("Logistic Coefficients");
    str(BetaLogit);
    
    if(i>1){
      cat("  LogLikComplete:",LogLikEvo[i], "Change:",LogLikEvo[i]-LogLikEvo[i-1]," Max Parameter Change:",max(MaxParamChange),"\n")
      if( (abs(LogLikEvo[i]-LogLikEvo[i-1])<0.02) && (max(MaxParamChange)<1E-3) ){
        ConvergedEM <- TRUE
      }
    }
    if((i>=sweeps)||ConvergedEM){
      IterateEM <- FALSE;
    }
  }
##############################################################################################

  ##browser()

  ## BY parts... 
  myNegBinIndRatios <- mapply(compNegBinLogRatio,  Rlist,NegBinParList)
  myMultinIndRatios <- mapply(compMultiNomLogRatio,Xlist,LambdaParList)
  myNegBinIndRatios <- diag(var(myNegBinIndRatios));
  names(myNegBinIndRatios) <- paste("Nb",names(myNegBinIndRatios))
  myMultinIndRatios <- diag(var(myMultinIndRatios));
  names(myMultinIndRatios) <- paste("Mn",names(myMultinIndRatios))
  myPriorLogRatios <- var(PriorLogRatio);
  names(myPriorLogRatios) <- "Prior"
  cc <- c(myPriorLogRatios,myNegBinIndRatios,myMultinIndRatios)
  cc <- cc/sum(cc);
  ##names(cc)[1] <- "Prior"
  cat("Model Parts Variance:\n");
  print(round(cc/sum(cc)*100,digits=2))

  myNegBinIndRatios <- mapply(compNegBinLogRatio,  Rlist,NegBinParList)
  myMultinIndRatios <- mapply(compMultiNomLogRatio,Xlist,LambdaParList)
  myNegBinIndRatios <- cor(LogRatios,myNegBinIndRatios);
  names(myNegBinIndRatios) <- paste("Nb",names(Rlist))
  myMultinIndRatios <- cor(LogRatios,myMultinIndRatios);
  names(myMultinIndRatios) <- paste("Mn",names(Rlist))
  myPriorLogRatios <- cor(LogRatios,PriorLogRatio);
  names(myPriorLogRatios) <- "Prior"
  dd <- c(myPriorLogRatios,myNegBinIndRatios,myMultinIndRatios)
  ##dd <- dd/sum(dd);
  cat("Model Parts Correlation:\n");
  print(round(dd*100,digits=2))
  
  ##browser()
  ##  list(PostPr=Ez,LogRatios=LogRatios,Lambda=Lambda, BetaLogit=BetaLogit,PriorPr=plogis(Y %*% BetaLogit),LogLikEvo=LogLikEvo, NBparams=c((A0*(1-B0))/B0, (A0*(1-B0))/(B0^2), (A1*(1-B1))/B1, (A1*(1-B1))/((B1)^2)),A0=A0,B0=B0,A1=A1,B1=B1) 
  list(PostPr=Ez,LogRatios=LogRatios,PriorLogRatio=PriorLogRatio,MultiNomLogRatio=MultiNomLogRatio,NegBinLogRatio=NegBinLogRatio
       ,LambdaParList=LambdaParList, BetaLogit=BetaLogit,PriorPr=plogis(Y %*% BetaLogit),LogLikEvo=LogLikEvo,LogLikEnd=LogLikEvo[i],NumIter=i
       ,NegBinParList=NegBinParList, ModelPartsVar=cc, ModelPartsCor=dd
       ,LogLikNoMix=LogLikNoMix,LogLikNoMixNoNegBin=LogLikNoMixNoNegBin
       ,NegBinParNoMixList=NegBinParNoMixList, LambdaParNoMixList=LambdaParNoMixList) 

}

