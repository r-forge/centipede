plotDistrAndCdf2 <-
function (Rpos,Rneg,main="Distribution of the reads",xlab="Number of reads",breaks){
  if(missing(breaks)){
    myMin=min(min(Rpos),min(Rneg))
    myMax=max(max(Rpos),max(Rneg))
    breaks=c(seq(myMin,myMax,(myMax-myMin+1)/50),1E9)
  }
  par(mfrow=c(2,1))
  par(mai=c(0.02,0.82,0.82,0.42))
  NegBars<-hist(Rneg,breaks=breaks, plot = FALSE)
  PosBars<-hist(Rpos,breaks=breaks, plot = FALSE)
##  barplot(rbind(PosBars$density,NegBars$density),names=format(NegBars$mids,digits=4),beside=T,ylab="Density",col=c("darkblue","lightblue"),legend=c("Positive","Negative"))
  barplot(rbind(PosBars$density/sum(PosBars$density)*100,NegBars$density/sum(NegBars$density)*100),names=format(NegBars$mids,digits=4),beside=T,ylab="Frequency %",col=c("darkblue","lightblue"),border=NA,legend=c("Positive","Negative"))
  title(main=main)
  NegEcdf<-ecdf(Rneg)
  PosEcdf<-ecdf(Rpos)
  par(mai=c(1.02,0.82,0.82,0.42))
  plot(PosEcdf,col.hor="darkblue",col.points="darkblue",cex=0.5,ylab="ECDF",xlab=xlab,main="");
  plot(NegEcdf,col.hor="lightblue",col.points="lightblue",add=TRUE,cex=0.5);
  ##ks.test(Rneg,Rpos,alternative="greater")
  par(mfrow=c(1,1))
}

