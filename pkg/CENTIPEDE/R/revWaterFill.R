revWaterFill <- function(Profile,cStop=0.25,xLim){
 center <- which.max(Profile)
 cMass=Profile[center]
 ii <- c(center-1,center+1)
 PP <- Profile[ii]
 cMass
 count <- 1;
 while(cMass<cStop){
##  print(PP)
  myMin <- which.max(PP)
##  print(myMin)
  if((myMin==1)){
    if(ii[1]<xLim[1]){
      PP[1] <- 0;
    }else{
      cMass <- cMass+PP[1]
      ii[1] <- ii[1]-1;
      PP[1] <- Profile[ii[1]]
    }
  }
  if((myMin==2)){
    if(ii[2]>xLim[2]){
      PP[2] <- 0;
    }else{
      cMass <- cMass+PP[2]
      ii[2] <- ii[2]+1;
      PP[2] <- Profile[ii[2]]
    }
  }
##
##  print(cMass)
##  print(ii)
  count <- count+1;
  if(count>length(Profile)*cStop)
    break;
  if(sum(PP)==0)
    break;
  }
 print(cMass)
 ii
}
