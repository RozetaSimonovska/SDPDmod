f_SIG<-function(rho,bhat,sige,W,Z,n,t,kz)
{
  Sn <- diag(n)-rho*W; Sni <- Matrix::solve(Matrix::Matrix(Sn, sparse = TRUE))
  Sni<-as.matrix(Sni)
  Gn <- W%*%Sni
  pterm <- sum(diag(Gn%*%Gn + Gn%*%t(Gn)))  #
  Gnt <- kronecker(diag(t),Gn)
  SIG <- matrix(0,nrow=(kz+2),ncol=(kz+2))   ###information matrix
  SIG[1:kz,1:kz] <- (1/(sige))*(t(Z)%*%Z)    # bhat,bhat
  SIG[1:kz,kz+1] <- (1/(sige))*(as.vector(t(Z)%*%Gnt%*%Z%*%t(t(bhat)))) #bhat,rho
  SIG[kz+1,1:kz] <- t(SIG[1:kz,kz+1])
  SIG[kz+1,kz+1] <- (1/(sige))*(as.vector(t(bhat)%*%t(Z)%*%t(
    as.matrix(Gnt))%*%as.matrix(Gnt)%*%Z%*%t(t(bhat)))) + pterm*t  #rho,rho
  SIG[kz+2,kz+2] <- (n*t)/(2*sige*sige)     #sige,sige
  SIG[kz+1,kz+2] <- (t/(sige))*sum(diag(Gn))
  SIG[kz+2,kz+1] <- SIG[kz+1,kz+2]  #rho,sige
  return(result=list(SIG=SIG,Gn=Gn,Sni=Sni))
}

f_OMG<-function(sige,Gn,n,kz,mu4){
  OMG <- matrix(0,nrow=(kz+2),ncol=(kz+2))
  OMG[kz+1,kz+1] <- c(diag(Gn)%*%diag(Gn)/n) ###
  OMG[kz+1,kz+2] <- (1/(2*sige*n))*sum(diag(Gn))  #
  OMG[kz+2,kz+1] <- OMG[kz+1,kz+2]
  OMG[kz+2,kz+2] <- 1/(4*sige*sige)
  OMG <- ((mu4-3*sige^2)/sige^2)*OMG
  return(OMG)
}
