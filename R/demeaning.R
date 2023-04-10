
demean<-function(y,x,N,t,effect,sind,tind){

  nobs<-N*t
  nvar<-ncol(x)

  meanny<-rep(0,N)
  meannx<-matrix(0,nrow=N,ncol=nvar)
  meanty<-rep(0,t)
  meantx<-matrix(0,nrow=t,ncol=nvar)

  if(effect %in% c("individual","twoways")){
    meanny<-tapply(y,sind,mean)
    for(i in 1:nvar){
      meannx[,i]<-tapply(x[,i],sind,mean)
    }
  }

  if(effect %in% c("time","twoways")){
    meanty<-tapply(y,tind,mean)
    for(i in 1:nvar)  meantx[,i]<-tapply(x[,i],tind,mean)
  }

  xtm<-matrix(0,nobs,nvar); xnm<-matrix(0,nobs,nvar); xmm<-matrix(0,nobs,nvar)
  for (i in 1:nvar)   xnm[,i]<-rep(meannx[,i],t)
  for (i in 1:nvar)   xtm[,i]<-rep(meantx[,i],each=N)
  for (i in 1:nvar)   xmm[,i]<-rep(mean(x[,i]),nobs)

  if(effect=="individual"){
    ywith<-y-rep(meanny,t)
    xwith<-x-xnm
  }else if(effect=="time"){
    ywith<-y-rep(meanty,each=N)
    xwith<-x-xtm
  }else if(effect=="twoways"){
    ywith<-y-rep(meanny,t)-rep(meanty,each=N)+rep(mean(y),nobs)
    xwith<-x-xnm-xtm+xmm
  }else if(effect=="none"){
    ywith<-y
    xwith<-x
  }else message("error effect")

  rez<-list(ywith,xwith,meanny,meannx,meanty,meantx)
  names(rez)<-c("yw","xw","mny","mnx","mty","mtx")
  return(rez)

}


demeanF<-function(y,x,N,t,effect,W){
  yf<-y
  xf<-x
  nt<-N*t; nv<-N; tv<-t
  k<-ncol(x)

  if(effect=="individual"){
    Jt<-Matrix::Matrix(diag(tv),
                       sparse = TRUE)-matrix(1/tv,
                                             nrow = tv,
                                             ncol = tv)
    Ftt<-eigen(Jt)$vectors
    f<-Ftt[,1:(tv-1)]
    yf<-matrix(yf, nrow = nv, ncol = tv)
    yf<-yf%*%f
    yf<-as.vector(yf)
    if(length(xf) != 0){
      xf<-array(xf,c(nv,tv,k))
      xtemp<-array(0,c(nv,tv-1,k))
      for(i in 1:k){
        xtemp[,,i]<-xf[,,i]%*%f
      }
      nt<-nt-nv
      xf<-matrix(xtemp, nrow = nt, ncol = k)
    } else {xf<-vector()}
    tv<-tv-1
    Wf<-W
  } else if(effect=="twoways"){
    Jt<-Matrix::Matrix(diag(tv),
                       sparse = TRUE)-matrix(1/tv,
                                             nrow = tv,
                                             ncol = tv)
    Ftt<-eigen(Jt)$vectors
    f<-Ftt[,1:(tv-1)]
    yf<-matrix(yf, nrow = nv, ncol = tv)
    yf<-yf%*%f
    if(length(xf)!=0){
      xf<-array(xf,c(nv,tv,k))
      xtemp<-array(0,c(nv,tv-1,k))
      for(i in 1:k){
        xtemp[,,i]<-xf[,,i]%*%f
      }
      nt<-nt-nv
      xf<-matrix(xtemp, nrow = nt, ncol = k)
    } else {xf<-vector()}
    tv<-tv-1

    Jn<-Matrix::Matrix(diag(nv),
                       sparse = TRUE)-matrix(1/nv,
                                             nrow = nv,
                                             ncol = nv)
    Fnn<-eigen(Jn)$vectors
    f<-Fnn[,1:(nv-1)]
    yf<-t(f)%*%yf
    yf<-matrix(yf, nrow = tv*(nv-1), ncol = 1)
    if(length(xf) != 0){
      xf<-array(xf,c(nv,tv,k))
      xtemp<-array(0,c(nv-1,tv,k))
      for(i in 1:k){
        xtemp[,,i]<-t(f)%*%xf[,,i]
      }
      nt<-nt-tv
      xf<-matrix(xtemp, nrow = nt, ncol = k)
    }else {xf<-vector()}
    nv<-nv-1
    Wf<-t(f)%*%W%*%f
  } else if(effect=="time"){
    yf<-y
    xf<-x
    nt<-N*t; nv<-N; tv<-t

    Jn<-Matrix::Matrix(diag(nv),
                       sparse = TRUE)-matrix(1/nv,
                                             nrow = nv,
                                             ncol = nv)
    Fnn<-eigen(Jn)$vectors
    f<-Fnn[,1:(nv-1)]
    yf<-matrix(yf, nrow = nv, ncol = tv)
    yf<-t(f)%*%yf
    yf<-matrix(yf, nrow = tv*(nv-1), ncol = 1)
    if(length(xf) != 0){
      xf<-array(xf,c(nv,tv,k))
      xtemp<-array(0,c(nv-1,tv,k))
      for(i in 1:k){
        xtemp[,,i]<-t(f)%*%xf[,,i]
      }
      nt<-nt-tv
      xf<-matrix(xtemp, nrow = nt, ncol = k)
    }else {xf<-vector()}

    nv<-nv-1
    Wf<-t(f)%*%W%*%f
  }
  rez<-list(yf,xf,Wf,nv,tv)
  names(rez)<-c("yf","xf","Wf","nv","tv")
  return(rez)
}


