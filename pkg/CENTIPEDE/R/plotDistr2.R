plotDistr2 <- function (Rpos,Rneg,main="Distribution of the reads",xlab="Number of reads",breaks,axes=TRUE,pos.legend="topright"){
  if(missing(breaks)){
    myMin=min(min(Rpos),min(Rneg))
    myMax=max(max(Rpos),max(Rneg))
    breaks=unique(c(seq(myMin,myMax,(myMax-myMin+1)/50),myMax))
  }
  par(cex=1.0)
  NegBars<-hist(Rneg,breaks=breaks, plot = FALSE)
  PosBars<-hist(Rpos,breaks=breaks, plot = FALSE)
  ##barplot(rbind(PosBars$density,NegBars$density),names=format(NegBars$mids,digits=4),beside=T,ylab="Density",col=c("darkblue","lightblue"),legend=c("Positive","Negative"))
  barplot(rbind(PosBars$counts/sum(PosBars$counts)*100,
                NegBars$counts/sum(NegBars$counts)*100),
          names=format(NegBars$mids,digits=4),axes=T,axisnames=F,
          beside=T,
          col=c("darkblue","lightblue"),
          border=NA,legend.text=c("Positive","Negative"),args.legend=list(x=pos.legend))
  if(axes)
    axis(1,at=seq(1,length(breaks)*3,3)-1,labels=round(breaks,digits=2))
  title(main=main,xlab=xlab,ylab="Frequency %",cex=1.5)
  ##lines(density(Rneg),col='lightblue',lwd=2)
  ##lines(density(Rpos),col='darkblue',lwd=2)
}
