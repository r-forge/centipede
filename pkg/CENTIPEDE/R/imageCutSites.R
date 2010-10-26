imageCutSites <- function (data,main="",xlab=""){
  data <- t(as.matrix(data));
  data[data>7] <- 7;
  data[data<0] <- 0;
  image(data,ylim=c(-0.00251,1.00251),col=c("white","darkblue","darkgreen","green","orange","yellow","red","darkred"),zlim=c(0,7))
  abline(v=c(0.25,0.75),lty=2,lwd=2)
  abline(v=c(0.5),lty=1,lwd=2)
  title(main=main,xlab=xlab)
  legend("bottomright",bg="white",c("0","1","2","3","4","5","6",">6"),fill=c("white","darkblue","darkgreen","green","orange","yellow","red","darkred"))
}

imageCutSitesCombined <- function (data,main="",xlab="",breaks){
  S <- dim(data)[2];
  data <- data[,1:(S/2)]+data[,(S/2+1):S];
  profile <- colMeans(data);
  data <- t(as.matrix(data));
  data[data<0] <- 0;
  ##data[data>7] <- 7;
  dd <- dim(data);
  if(missing(breaks)){
    #breaks=breaks=c(0,1,2,4,8,16,32,64,1E9)
    breaks=c(0:7,1E9)
  }
  data <- cut(data,breaks=breaks,right=F);
  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,.2)))(length(breaks)-1);  
  myNames <- levels(data);
  data <- matrix(as.numeric(data)-1,dd)
  image(data,ylim=c(-0.00251,1.00251),col=myColors,zlim=c(0,7),axes=FALSE)
  box()
  abline(v=0.5,lty=2,lwd=1)
##  lines((1:length(profile))/length(profile),profile*0.5/max(profile),lwd=2,col='red')
  title(main=main,xlab=xlab)
##  legend("right",bg="white",myNames,fill=myColors)
}
