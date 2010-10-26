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
fitNegBinMixture <- function(R,K,p,NegBinParList,sweeps=50,TrimP=0.0001,NRiter=5) {
  LogLikEvo <- vector(length=sweeps)

  ## Checking Number of Matrices... and dimensionality match
  
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


  ## Initialization of the Negative Binomial, looks a bit horrible...
  ## But it is simply method of moments plus some clipping to aboid extreme values. 
  initNegBinParams <- function(R,Ez,DampNegBin=0.0){
    Ez <- Ez*(1-DampNegBin)+0.5*(DampNegBin);
    myM<-sum(R*Ez)/sum(Ez)
    myV<-sum(R^2 * Ez)/sum(Ez) - myM^2
    B1<-myM/myV;
        if(B1 > 0.999999){
        B1 <- 0.999999   
        }      
    A1<-myM*B1/(1-B1);
    list(logA=log(A1),logitB=qlogis(B1))
  }

  compNegBinParamsMixture <- function(R,NegBinPar,Ez,Pcut=1,NRiter=5){
    ## Negative Binomial dampening experiment!
    ## Fast update first
      NegBinPar$logitB <- NegBinPar$logA + log(sum(Ez)) - log(sum(R*Ez))    
      NegBinPar1 <- c(NegBinPar$logA,NegBinPar$logitB)
      ## M-Step NegBinom parameters 
      NegBinPar1 <- optim(NegBinPar1, fNegBin2, gNegBin2, R=R, Ez=Ez , Pcut=Pcut, method="L-BFGS-B",lower=c(-20,-20),upper=c(20,20),control=list(maxit=NRiter))    
      NegBinPar$logA <- NegBinPar1$par[1];
      NegBinPar$logitB <- NegBinPar1$par[2];      
      NegBinPar
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


  compNegBinLogLik <- function(R,NegBinPar,Pcut=1,DampNegBin=0.0){
    A <- exp(NegBinPar$logA)
    logitB <- NegBinPar$logitB-log(Pcut)
    NegBinLogRatio <-  lgamma(R+A) - lgamma(A) + A*logitB - (R+A)*log(1+exp(logitB));
    NegBinLogRatio
  }
  

################################################################
################################################################
## Initializiation / or checking given Intialization..
################################################################  
  cat("Initialization of the parameters:\n")

  L <- length(R)
  
  if(missing(K)){
    K <- length(p)
    stopifnot(length(NegBinParList)==K)
    stopifnot(sum(p)==1)
  }else{    
    ## NegBin initialization...   
    ## M-Step NegBinParameters
    
    ## assert missing p, and NegBinParList...
    mybreaks=quantile(R,seq(0,1,len=2+K-1))
    Pmat <- matrix(0,K,L);
    aa <- cut(R,mybreaks,include.lowest=TRUE,labels=F)
    Pmat[((1:L) - 1) * K + aa] <- 1    
    Pmat <- Pmat*0.9+0.1*1/K
    Pmat <- t(Pmat)
    aa <- hist(R,breaks=mybreaks,plot=F)$counts
    p <- aa/sum(aa);
    NegBinParList <- lapply(1:K,function(kk){
      initNegBinParams(R,Ez=Pmat[,kk])})
    NegBinParList <- lapply(1:K,function(kk){
      compNegBinParamsMixture(R,NegBinParList[[kk]],Ez=Pmat[,kk])})    
  }
   
############################################################
  
 ## Since I will do the difference with the backgraound only model I don't need this anymore!...
 MyLogLikConst <- -sum(lgamma(R)) #( -log(S)*sum(X)); # - sum(lgamma(1+X)) - log(S)*sum(X)  #These guys are constants, so if we take them out, we go faster...

#########################################################################################################
##  EM loop starts here ##
#########################################################################################################
  ## First E-step outside the loop
  cat("\n Starting EM loop:\n")
    
  ##NegBinLogRatio <- mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin));
  NegBinLogLik <- sapply(1:K,function(kk){
    compNegBinLogLik(R,NegBinParList[[kk]])+log(p[kk])})
  ## NegBinLogRatio <- sapply(1:K,function(kk){
  ##   -log(rowSums(exp(NegBinLogLik-NegBinLogLik[,kk])))
  ## })
  Pmat <- sapply(1:K,function(kk){
    1/(rowSums(exp(NegBinLogLik-NegBinLogLik[,kk])))
  }) 
 
  LogLikEvo <- 0;

  aa <- rowMax(NegBinLogLik)
  LogLikEvo <- sum(log(rowSums(exp(NegBinLogLik-aa)))+aa)
  
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
    
    
    ## M-Step NegBin
    OldNegBinPars <- NegBinParList;
    NegBinParList <- lapply(1:K,function(kk){
      compNegBinParamsMixture(R,NegBinParList[[kk]],Ez=Pmat[,kk])})

    ## M-step Prior
    OldPrior <- p
    p <- colSums(Pmat)/L
    
    ##NegBinLogRatio <- mapply(compNegBinLogRatio,Rlist,NegBinParList,MoreArgs=list(TrimP=TrimP,DampNegBin=DampNegBin));
    NegBinLogLik <- sapply(1:K,function(kk){
      compNegBinLogLik(R,NegBinParList[[kk]])+log(p[kk])})
    Pmat <- sapply(1:K,function(kk){
      1/(rowSums(exp(NegBinLogLik-NegBinLogLik[,kk])))
    })

    aa <- rowMax(NegBinLogLik)
    LogLikEvo[i] <- sum(log(rowSums(exp(NegBinLogLik-aa)))+aa)   

    ## Parameter Convergence monitoring..... 
    NegBinChange <-as.matrix(mapply(function(i,j){(as.data.frame(j)-as.data.frame(i))},OldNegBinPars,NegBinParList))
    PriorChange <- (p-OldPrior)
    MaxParamChange <- c(max(abs(unlist(NegBinChange))),max(abs(PriorChange)))
    str(MaxParamChange)

    print("Neg Bin Changes:");
    print(NegBinChange);
    print("Neg Bin New Values:");
    print(sapply(NegBinParList,cbind))
    print("Prior Changes:");
    print(PriorChange);
    print("Prior Coefficients");
    ##str(BetaLogit);
    print(p)
    
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
 

  myR <- 0:max(R)
  NegBinLogLik <- sapply(1:K,function(kk){
    compNegBinLogLik(myR,NegBinParList[[kk]])+log(p[kk])})
  aa <- rowMax(NegBinLogLik)
  aux <- lgamma(myR+1)
  myP <- exp(log(rowSums(exp(NegBinLogLik-aa)))+aa-aux)
  ##
  emP <- hist(R,breaks=c(myR-0.5,1E20),plot=F)$counts
  emP <- emP/sum(emP)
  par(mfrow=c(2,2))
  ##plot(myR,emP,t='h',xlim=c(0,1000))
  plot(myR,emP,t='h',xlim=c(0,quantile(R,0.75)))
  mycol <- colorRampPalette(c('green','darkgreen'))(K)
  sapply(1:K,function(kk){lines(myR,exp(NegBinLogLik[,kk]-aux),col=mycol[kk])})
  lines(myR,myP,t='l',col='red',lwd=3)
  ##
  plot(myR,emP,t='h',log='x')
  mycol <- colorRampPalette(c('green','darkgreen'))(K)
  sapply(1:K,function(kk){lines(myR,exp(NegBinLogLik[,kk]-aux),col=mycol[kk])})
  lines(myR,myP,t='l',col='red',lwd=3)
  ##
  qx <- quantile(R,seq(0,1,len=200))+1
  px <- cumsum(myP)[qx]
  qy <- quantile(R,px)
  plot(qx,qy)
  abline(0,1,col='red')
  ##
  plot(ecdf(R))
  lines(myR,cumsum(myP),col='red',lwd=1)
  
  browser()
  
 ## list(PostPr=PiMat,LogRatios=LogRatios,PriorLogRatio=PriorLogRatio,MultiNomLogRatio=MultiNomLogRatio,NegBinLogRatio=NegBinLogRatio
 ##      ,Lambda=Lambda, LambdaUnc=LambdaUnc, BetaLogit=BetaLogit, PriorPr=PriorMat, LogLikEvo=LogLikEvo,LogLikEnd=LogLikEvo[i],NumIter=i
 ##      ,NegBinParList=NegBinParList) 
}

