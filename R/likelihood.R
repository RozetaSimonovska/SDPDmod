
f_sar<-function(rho,detval,epeo,eped,epeod,n,t,dynamic){
  gsize <- detval[2,1]-detval[1,1]
  i1 <- max(which(detval[,1] <= rho + gsize))
  i2 <- max(which(detval[,1]<= rho - gsize))
  ind <- round((i1+i2)/2)
  if(length(ind)==0){  ind <- 1}
  detm <- detval[ind,2]
  zz <- epeo - 2*rho*epeod + rho*rho*eped
  if(dynamic){  zz <- zz/(n*t)}
  llike <- (n*t/2)*log(zz) - t*detm
  return(llike)
}

f2_sar<-function(rho,bhat,y,Z,Wy,detval,n,t,sige){
  gsize <- detval[2,1] - detval[1,1]
  i1 <- max(which(detval[,1] <= rho + gsize))
  i2 <- max(which(detval[,1]<= rho - gsize))
  indr <- round((i1+i2)/2)
  if(length(indr)==0){  indr <- 1}
  detm <- detval[indr,2]
  eee <- y-Z%*%bhat-rho*Wy
  epee <- c(t(eee)%*%eee)
  llik<- as.numeric((-(n*t/2)*log(2*pi*sige)) + t*detm -( 1/(2*sige))*epee)
  return(llik)
}

f2_sar2<-function(rho,bhat,y0,Z0,detval,n0,t0,sige,sind0,tind0,Wy0){
  x<-cbind(Wy0,Z0)
  re1<-demean(y0,x=x,N=n0,t=t0,effect="twoways",sind=sind0,tind=tind0)
  y<-re1$yw; Xx<-re1$xw
  n<-n0; t<-t0
  Wy<-Xx[,1]; if(length(Z0)!=0){ Z<-Xx[,2:ncol(x)] }

  gsize <- detval[2,1] - detval[1,1]
  i1 <- max(which(detval[,1] <= rho + gsize))
  i2 <- max(which(detval[,1]<= rho - gsize))
  indr <- round((i1+i2)/2)
  if(length(indr)==0){  indr <- 1}
  detm <- detval[indr,2]
  eee <- y-Z%*%bhat-rho*Wy
  epee <- c(t(eee)%*%eee)
  llik<- as.numeric((-(n*t/2)*log(2*pi*sige)) + t*detm -( 1/(2*sige))*epee)
  return(llik)
}

f2_sar_dyn<-function(rho,bhat,y0,Z0,detval,n0,t0,sige,sind0,tind0,Wy0){
  x<-cbind(Wy0,Z0)
  re1<-demean(y0,x=x,N=n0,t=t0,effect="individual",sind=sind0,tind=tind0)
  y<-re1$yw; Xx<-re1$xw
  n<-n0; t<-t0
  Wy<-Xx[,1]; if(length(Z0) != 0) Z<-Xx[,2:ncol(x)]

  gsize <- detval[2,1] - detval[1,1]
  i1 <- max(which(detval[,1] <= rho + gsize))
  i2 <- max(which(detval[,1]<= rho - gsize))
  indr <- round((i1+i2)/2)
  if(length(indr)==0){  indr <- 1}
  detm <- detval[indr,2]
  eee <- y-Z%*%bhat-rho*Wy
  epee <- c(t(eee)%*%eee)
  llik<- as.numeric((-(n*t/2)*log(2*pi*sige)) + t*detm -( 1/(2*sige))*epee)
  return(llik)
}
