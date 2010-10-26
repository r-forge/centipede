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
fitMultiCentipede <- function(Xlist,Y, KmerSeq,Mappability, KmerCutPr, BetaLogit, Lambda, NegBinParList,sweeps=50,DampLambda=0.01,DampNegBin=0.0,TrimP=0.0001,NRiter=5,BetaLogitLower=-40,BetaLogitUpper=20,FixLambda=FALSE) {
  LogLikEvo <- vector(length=sweeps)

  ## Checking Number of Matrices... and dimensionality match
  stopifnot(is.list(Xlist)); 
  NumTissues <- length(Xlist); 
  stopifnot(NumTissues>0);
   
  ## Y<-as.matrix(Y)
  MatrixMode <- prod(sapply(Xlist,is.matrix))==1
  BigMatrixMode <- prod(sapply(Xlist,is.big.matrix))==1
  stopifnot(MatrixMode || BigMatrixMode);
      
  L <- dim(Xlist[[1]])[1]  #Number of rows has to be the same on all matrices
  stopifnot(sapply(Xlist,function(i) {dim(i)[1]==L}))

  #Number of columns cannot differ for the multi-tissue model
  S <- dim(Xlist[[1]])[2]  #Number of cols has to be the same on all matrices
  stopifnot(sapply(Xlist,function(i) {dim(i)[2]==S}))

  if(missing(Y)){
    ##Y <- matrix(1,L); ## Intercept point.
    Y <- matrix(0,L,0);
  }
  stopifnot(is.matrix(Y));
  stopifnot(dim(Y)[1]==L);
  
  NumCov<-dim(Y)[2]

  if(MatrixMode)
    Rlist <- lapply(Xlist,rowSums) #Total number of reads
  if(BigMatrixMode)
    Rlist <- lapply(Xlist,function(X){rowSums(as.matrix(X))}) #Total number of reads
    ##Rlist <- lapply(Xlist,function(X){apply(X,1,sum)}) #Total number of reads    
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
  compPriorLogRatio <- function(Y,BetaLogit){
    C <- dim(Y)[2]
    K <- length(BetaLogit)-C;
    BetaInt <- BetaLogit[1:K]
    BetaCom <- BetaLogit[(1:C)+K]
    Cl <- Y %*% BetaCom;
    Clk <- sapply(1:K,function(kk){
      Cl + BetaInt[kk]});
    Clk
  }

  fLogitMulti <- function(BetaLogit,Y,PiMat){
    C <- dim(Y)[2]
    K <- dim(PiMat)[2]
    BetaInt <- BetaLogit[1:K]
    Cl <- 0;
    if(C>0){
      BetaCom <- BetaLogit[(1:C)+K]
      Cl <- Y %*% BetaCom;
    }
    fk <- sapply(1:K,function(kk)
                 {
                   Clk <- Cl + BetaInt[kk];
                   sum(PiMat[,kk] * Clk) - sum(log(1 + exp(Clk)));
                 });
    -sum(fk);
  }
  gLogitMulti <- function(BetaLogit,Y,PiMat){
    C <- dim(Y)[2]
    K <- dim(PiMat)[2]
    BetaInt <- BetaLogit[1:K]
    Cl <- 0;
    if(C>0){
      BetaCom <- BetaLogit[(1:C)+K]
      Cl <- Y %*% BetaCom;
    }
    Clk <- sapply(1:K,function(kk){
      Cl + BetaInt[kk]});
    Plk <- plogis(Clk);
    BetaIntGr <- colSums(PiMat-Plk)
    BetaComGr <- colSums(t(PiMat-Plk) %*% Y)
    -c(BetaIntGr,BetaComGr)
  }

  
  ## Preparing functions for the Negative Binomial alpha_0 using R's BGFS -- 
  ## Objective function
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

  ## Objective function
  fMultSeqDepMultiTissue <- function(LogPar,Xsk,EzRlk,rho){
    lambda <- exp(LogPar);
    lambda <- lambda/sum(lambda);
    ##aux <- apply(rho,1,function(rhorow){sum(rhorow*lambda)});
    aux <- rho %*% lambda
    -(sum(log(lambda) %*% Xsk) - sum(t(log(aux)) %*% EzRlk))
  }
  gMultSeqDepMultiTissue <- function(LogPar,Xsk,EzRlk,rho){
    lambda <- exp(LogPar);
    lambda <- lambda/sum(lambda);
    ##aux <- apply(rho,1,function(rhorow){sum(rhorow*lambda)});
    aux <- rho %*% lambda
    ##-(Xs-lambda*colSums(EzRl/aux*rho))
    -rowSums(sapply(1:dim(Xsk)[2],function(kk){(Xsk[,kk]-lambda*t(t(EzRlk[,kk]/aux) %*% rho))}))
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
  
  compLambdaMultiTissue <- function(Xlist,PiMat,DampLambda=0.0){
    K <- length(Xlist)
    if(MatrixMode)
      Lambda <- rowSums(sapply(1:K,function(kk){(t(PiMat[,kk]) %*%  Xlist[[kk]])}))
    if(BigMatrixMode)
      Lambda <- rowSums(sapply(1:K,function(kk){(t(PiMat[,kk]) %*%  as.matrix(Xlist[[kk]]))}))
    Lambda <- (Lambda)/sum(Lambda);
    Lambda <- Lambda*(1-DampLambda) + (DampLambda)*1/S;
    Lambda
  }

  compLambdaSeqDepMultiTissue <- function(LambdaIni,Xlist,Rlist,PiMat,DampLambda=0.0,rho){
    K <- length(Xlist)
    LogParIni <- log(LambdaIni)
    ##Xsk,EzRlk
    if(MatrixMode)
      Xsk = sapply(1:K,function(kk){t(t(PiMat[,kk]) %*% Xlist[[kk]])})
    if(BigMatrixMode)
      Xsk = sapply(1:K,function(kk){t(t(PiMat[,kk]) %*% as.matrix(Xlist[[kk]]))})    
    EzRlk = sapply(1:K,function(kk){PiMat[,kk] * Rlist[[kk]]})
    optPar <- optim(LogParIni,fMultSeqDepMultiTissue,gMultSeqDepMultiTissue,
                      Xsk = Xsk, EzRlk = EzRlk, rho = rho,
                      method="BFGS",control = list(maxit = NRiter))
##    str(optPar)
    LogPar <- optPar$par
    Lambda <- exp(LogPar);
    ##sum(Lambda)
    Lambda <- Lambda/sum(Lambda);
    Lambda <- Lambda*(1-DampLambda) + (DampLambda)*1/S;
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
      if(MatrixMode)              
        MultiNomLogRatio <-  X %*% log(Lambda*S); ## Multinomial for 1 and 0
      if(BigMatrixMode)              
        MultiNomLogRatio <-  as.matrix(X) %*% log(Lambda*S); ## Multinomial for 1 and 0      
      ##MultiNomLogRatio <- compMultLogRatioOccPmax2(X,Lambda)      
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
      if(MatrixMode)
        MultiNomLogRatio <-  X %*% log(Lambda) - R*logaux; ## Multinomial for 1 and 0
      if(BigMatrixMode)
        MultiNomLogRatio <-  as.matrix(X) %*% log(Lambda) - R*logaux; ## Multinomial for 1 and 0      
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
    if(MatrixMode)
      aux <- rowSums(X * log(rho));
    if(BigMatrixMode)
      aux <- rowSums(as.matrix(X) * log(rho));    
    if(DampNegBin<1){
      A0 <- exp(NegBinPar$logA0)
      logitB0 <- NegBinPar$logitB0
      LogLik0 <-  lgamma(R+A0) - lgamma(A0) + A0*logitB0 - (R+A0)*log(1+exp(logitB0)) + aux;
      LogLik0 <- clipExtremes(LogLik0,TrimP)
    }else{
      LogLik0 <- aux;
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
## Initializiation / or checking given Intialization..
################################################################  
  cat("Initialization of the parameters:\n")

  ## Could add all-tissues together...
  PiMat <- sapply(1:NumTissues,function(k){(Rlist[[k]]> quantile(Rlist[[k]],0.9))+0.0})
  
  ## BetaLogit initialization, ....  
  if(missing(BetaLogit)){
    BetaLogit <- matrix(0,NumTissues+NumCov,1)
    BetaLogitLower <- rep(-25,NumTissues+NumCov)
    BetaLogitLower[NumTissues+1] <- 1;
    BetaLogitUpper <- rep(25,NumTissues+NumCov)
    BetaLogitUpper[1:NumTissues] <- -3;
    BetaFit <- optim(BetaLogit,fLogitMulti,gLogitMulti,Y=Y,PiMat=PiMat,
                     method="L-BFGS-B",lower=BetaLogitLower,upper=BetaLogitUpper);
    BetaLogit<-BetaFit$par;    
  }
  
  row.names(BetaLogit) <- c(names(Xlist),colnames(Y));

  #PriorLogRatio <- Y %*% BetaLogit;
  #Ez<-plogis(PriorLogRatio);

  ## Lambda initialization...
  ## M-Step Lambda 
  if(missing(Lambda)){
    stopifnot(FixLambda==FALSE);
    Lambda <- compLambdaMultiTissue(Xlist,PiMat,DampLambda=DampLambda)
  ##  LambdaParList <- lapply(Xlist,compLambda,Ez=Ez,DampLambda=DampLambda);
  }
  
  if(SeqDepFlag){
    Lambda <- compLambdaSeqDepMultiTissue(Lambda,Xlist,Rlist,PiMat,DampLambda=DampLambda,rho=RhoMat);
  }
  
  ## NegBin initialization...   
  ## M-Step NegBinParameters
  if(missing(NegBinParList)){
    NegBinParList <- lapply(1:NumTissues,function(kk){
      initNegBinParamsMixture(Rlist[[kk]],Ez=PiMat[,kk],DampNegBin=DampNegBin)})
    NegBinParList <- lapply(1:NumTissues,function(kk){
      compNegBinParamsMixture(Rlist[[kk]],NegBinParList[[kk]],Ez=PiMat[,kk],DampNegBin=DampNegBin)[[1]]})
  }
   names(NegBinParList) <- names(Xlist)
  
  ##########################
  ## No mixture Model Fit ##
  ##########################
  
##  NegBinParNoMixList <- lapply(Rlist,compNegBinParams,NRiter=100)
##  LambdaParNoMixList <- lapply(Xlist,compLambdaNoMix,DampLambda=DampLambda)
##  LogLikNoMix <- sum(mapply(compLogLik1,Xlist,LambdaParNoMixList,NegBinParNoMixList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)))
##  LogLikNoMixNoNegBin <- sum(mapply(compLogLik1,Rlist,LambdaParNoMixList,NegBinParNoMixList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)))
  
############################################################
  
 ## Since I will do the difference with the backgraound only model I don't need this anymore!...
 MyLogLikConst <- 0 #( -log(S)*sum(X)); # - sum(lgamma(1+X)) - log(S)*sum(X)  #These guys are constants, so if we take them out, we go faster...

#########################################################################################################
##  EM loop starts here ##
#########################################################################################################
  ## First E-step outside the loop
  cat("\n Starting EM loop:\n")
    
  NegBinLogRatio <- mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin));
  if(SeqDepFlag){
    MultiNomLogRatio <- mapply(compMultiNomLogRatioSeqDep,Xlist,Rlist,MoreArgs=list(Lambda=Lambda,rho=RhoMat, TrimP=TrimP));
  }else{
    MultiNomLogRatio <- mapply(compMultiNomLogRatio,Xlist,MoreArgs=list(Lambda=Lambda,TrimP=TrimP));
  }
  ##PriorLogRatio <- Y %*% BetaLogit;
  PriorLogRatio <- compPriorLogRatio(Y,BetaLogit)
  ## Plk <- plogis(Clk);

  LogRatios <- PriorLogRatio+NegBinLogRatio+MultiNomLogRatio;
  PiMat<-plogis(LogRatios);

  LogLikEvo <- 0;
  if(SeqDepFlag){
    LogLik0 <- sum(mapply(compLogLik0SeqDep,Xlist,Rlist,NegBinParList,MoreArgs=list(rho=RhoMat,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
  }else{
    LogLik0 <- sum(mapply(compLogLik0,Rlist,NegBinParList,MoreArgs=list(S=S,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
  }      
  Aux <- log(1+exp(LogRatios));
  ## LogLik0 <- sum(mapply(compLogLik0,Rlist,NegBinParList,MoreArgs=list(S=S,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));  
  Aux <- log(1+exp(LogRatios));
  Aux[LogRatios>18] <- LogRatios[LogRatios>18]; ## makes the thing more stable if LogRatio very large! so exp(LogRatio) does not explode
  LogLikEvo <- sum(Aux) + LogLik0; ## + sum(log(1-myPi)) - L*lgamma(A0) + L*log(B0)*A0 + sum(lgamma(R+A0)) + log(1-B0)*sum(R) + MyLogLikConst 
##  - sum(R*log(S))

#########################################################################################################
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
    
    ## M-Step Lambda    
    OldLambdas <- Lambda; ##LambdaParList;
    if(FixLambda==FALSE){
      if((i<10)|((i%%5)==0)){
        if(SeqDepFlag){
          Lambda <- compLambdaSeqDepMultiTissue(Lambda,Xlist,Rlist,PiMat,DampLambda=DampLambda,rho=RhoMat);      
        }else{
          Lambda <- compLambdaMultiTissue(Xlist,PiMat,DampLambda=DampLambda);
        }
      }
    }
    
    ## M-Step NegBin
    OldNegBinPars <- NegBinParList;
    ## ##NegBinParList <- lapply(Rlist,initNegBinParamsMixture,Ez=Ez);
    ##NegBinParList <- mapply(compNegBinParamsMixture,Rlist,NegBinParList,MoreArgs=list(Ez=Ez,NRiter=NRiter,DampNegBin=DampNegBin,Pcut=Pcut));
    NegBinParList <- lapply(1:NumTissues,function(kk){
      compNegBinParamsMixture(Rlist[[kk]],NegBinParList[[kk]],Ez=PiMat[,kk],DampNegBin=DampNegBin)[[1]]})


    ## M-step Logistic,
    OldBetaLogit <- BetaLogit;
##    BetaFit <- optim(BetaLogit,fLogit,gLogit,Y=Y,Ez=Ez,method="L-BFGS-B",lower=rep(-20,K),upper=rep(20,K),control=list(maxit=5));
##    BetaFit <- optim(BetaLogit,fLogit,gLogit,Y=Y,Ez=Ez,method="BFGS",control=list(maxit=NRiter));
##    BetaLogit<-BetaFit$par;
    ##BetaFit <- optim(BetaLogit,fLogitMulti,gLogitMulti,Y=Y,PiMat=PiMat,method="BFGS");
    BetaFit <- optim(BetaLogit,fLogitMulti,gLogitMulti,Y=Y,PiMat=PiMat,
                     method="L-BFGS-B",lower=BetaLogitLower,upper=BetaLogitUpper,control=list(maxit=NRiter));
    BetaLogit<-BetaFit$par;    


    ##myPi <- plogis(Y %*% BetaLogit)
    ##PriorLogRatio <- Y %*% BetaLogit;
    ## if(SeqDepFlag){
    ##   NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,Pcut=Pcut,DampNegBin=DampNegBin)));
    ##   ##browser();
    ##   MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatioSeqDep,Xlist,Rlist,LambdaParList,MoreArgs=list(rho=RhoMat, TrimP=TrimP)));
    ## }else{
    ##   NegBinLogRatio <- rowSums(mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin)));
    ##   MultiNomLogRatio <- rowSums(mapply(compMultiNomLogRatio,Xlist,LambdaParList,MoreArgs=list(TrimP=TrimP)));
    ## }
    ## LogRatios <- PriorLogRatio+NegBinLogRatio+MultiNomLogRatio;
    ## Ez<-plogis(LogRatios);

    NegBinLogRatio <- mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin));
    if(FixLambda==FALSE){
      if((i<10)|((i%%5)==0)){        
        if(SeqDepFlag){
          MultiNomLogRatio <- mapply(compMultiNomLogRatioSeqDep,Xlist,Rlist,MoreArgs=list(Lambda=Lambda,rho=RhoMat, TrimP=TrimP));
        }else{
          MultiNomLogRatio <- mapply(compMultiNomLogRatio,Xlist,MoreArgs=list(Lambda=Lambda,TrimP=TrimP));
        }
      }
    }
    ##PriorLogRatio <- Y %*% BetaLogit;
    PriorLogRatio <- compPriorLogRatio(Y,BetaLogit)
    ## Plk <- plogis(Clk);
    LogRatios <- PriorLogRatio+NegBinLogRatio+MultiNomLogRatio;
    PiMat<-plogis(LogRatios);

    if(SeqDepFlag){
      if((i<10)|((i%%5)==0)){
        LogLik0 <- sum(mapply(compLogLik0SeqDep,Xlist,Rlist,NegBinParList,MoreArgs=list(rho=RhoMat,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
        LogLik0Adj <- LogLik0-sum(mapply(compLogLik0,Rlist,NegBinParList,MoreArgs=list(S=S,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
      }else{
        LogLik0 <- LogLik0Adj+sum(mapply(compLogLik0,Rlist,NegBinParList,MoreArgs=list(S=S,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
      }
    }else{
      LogLik0 <- sum(mapply(compLogLik0,Rlist,NegBinParList,MoreArgs=list(S=S,TrimP=TrimP,DampNegBin=DampNegBin))) - sum(log(1+exp(PriorLogRatio)));
    }      
    Aux <- log(1+exp(LogRatios));
    Aux[LogRatios>18] <- LogRatios[LogRatios>18]; ## makes the thing more stable if LogRatio very large! so exp(LogRatio) does not explode
    LogLikEvo[i] <- sum(Aux) + LogLik0; ## + sum(log(1-myPi)) - L*lgamma(A0) + L*log(B0)*A0 + sum(lgamma(R+A0)) + log(1-B0)*sum(R) + MyLogLikConst 

    ## Parameter Convergence monitoring..... 
    ##LambdaChange <-sapply(mapply(function(i,j){abs(j-i)},OldLambdas,LambdaParList),max);
    LambdaChange <- max(abs(OldLambdas-Lambda))
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
    ##print("Neg Bin Changes:");
    ##print(NegBinChange);
    print("Neg Bin New Values:");
    print(sapply(NegBinParList,cbind))
    ##print("Logistic Changes:");
    ##print(LogitChange);
    print("Logistic Coefficients");
    ##str(BetaLogit);
    print(t(BetaLogit))
    
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

  

  PriorMat<-plogis(PriorLogRatio);
  LambdaUnc <- compLambdaMultiTissue(Xlist,PiMat,DampLambda=DampLambda)
  ## plot(Lambda,t='l')
  ## lines(LambdaUnc,t='l',col='red')
  ## plot(PiMat[,1] ,PiMat[,2],cex=0.3)
  ## points(PriorMat[,1] ,PriorMat[,2],cex=0.3,col='red')
  ## plot(PiMat[,1] ,PiMat[,4],cex=0.3)
  ## points(PriorMat[,1] ,PriorMat[,4],cex=0.3,col='red')
  ## browser()

  ## BY parts... 
##  myNegBinIndRatios <- mapply(compNegBinLogRatio,  Rlist,NegBinParList)
##  myMultinIndRatios <- mapply(compMultiNomLogRatio,Xlist,LambdaParList)
##  myNegBinIndRatios <- diag(var(myNegBinIndRatios));
##  names(myNegBinIndRatios) <- paste("Nb",names(myNegBinIndRatios))
##  myMultinIndRatios <- diag(var(myMultinIndRatios));
##  names(myMultinIndRatios) <- paste("Mn",names(myMultinIndRatios))
##  myPriorLogRatios <- var(PriorLogRatio);
##  names(myPriorLogRatios) <- "Prior"
##  cc <- c(myPriorLogRatios,myNegBinIndRatios,myMultinIndRatios)
##  cc <- cc/sum(cc);
  ## ##names(cc)[1] <- "Prior"
  ## cat("Model Parts Variance:\n");
  ## print(round(cc/sum(cc)*100,digits=2))

  ## myNegBinIndRatios <- mapply(compNegBinLogRatio,  Rlist,NegBinParList)
  ## myMultinIndRatios <- mapply(compMultiNomLogRatio,Xlist,LambdaParList)
  ## myNegBinIndRatios <- cor(LogRatios,myNegBinIndRatios);
  ## names(myNegBinIndRatios) <- paste("Nb",names(Rlist))
  ## myMultinIndRatios <- cor(LogRatios,myMultinIndRatios);
  ## names(myMultinIndRatios) <- paste("Mn",names(Rlist))
  ## myPriorLogRatios <- cor(LogRatios,PriorLogRatio);
  ## names(myPriorLogRatios) <- "Prior"
  ## dd <- c(myPriorLogRatios,myNegBinIndRatios,myMultinIndRatios)
  ## ##dd <- dd/sum(dd);
  ## cat("Model Parts Correlation:\n");
  ## print(round(dd*100,digits=2))
  
  ##browser()
  ##  list(PostPr=Ez,LogRatios=LogRatios,Lambda=Lambda, BetaLogit=BetaLogit,PriorPr=plogis(Y %*% BetaLogit),LogLikEvo=LogLikEvo, NBparams=c((A0*(1-B0))/B0, (A0*(1-B0))/(B0^2), (A1*(1-B1))/B1, (A1*(1-B1))/((B1)^2)),A0=A0,B0=B0,A1=A1,B1=B1) 
  ## list(PostPr=Ez,LogRatios=LogRatios,PriorLogRatio=PriorLogRatio,MultiNomLogRatio=MultiNomLogRatio,NegBinLogRatio=NegBinLogRatio
  ##      ,LambdaParList=LambdaParList, BetaLogit=BetaLogit,PriorPr=plogis(Y %*% BetaLogit),LogLikEvo=LogLikEvo,LogLikEnd=LogLikEvo[i],NumIter=i
  ##      ,NegBinParList=NegBinParList, ModelPartsVar=cc, ModelPartsCor=dd
  ##      ,LogLikNoMix=LogLikNoMix,LogLikNoMixNoNegBin=LogLikNoMixNoNegBin
  ##      ,NegBinParNoMixList=NegBinParNoMixList, LambdaParNoMixList=LambdaParNoMixList)
  list(PostPr=PiMat,LogRatios=LogRatios,PriorLogRatio=PriorLogRatio,MultiNomLogRatio=MultiNomLogRatio,NegBinLogRatio=NegBinLogRatio
       ,Lambda=Lambda, LambdaUnc=LambdaUnc, BetaLogit=BetaLogit, PriorPr=PriorMat, LogLikEvo=LogLikEvo,LogLikEnd=LogLikEvo[i],NumIter=i
       ,NegBinParList=NegBinParList) 
}

