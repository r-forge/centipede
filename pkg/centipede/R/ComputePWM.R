ComputePWM <-
function(Sequence,Ez){
  if(missing(Ez)){
    Ez <- rep(1,dim(Sequence)[1])
  }
  Am <- (t(Ez)%*%((Sequence=='A')+0.0));
  Cm <- (t(Ez)%*%((Sequence=='C')+0.0));
  Gm <- (t(Ez)%*%((Sequence=='G')+0.0));
  Tm <- (t(Ez)%*%((Sequence=='T')+0.0));         
  rbind(Am,Cm,Gm,Tm)/sum(Ez)
}

