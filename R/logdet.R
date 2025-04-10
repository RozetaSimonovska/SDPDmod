lndetfull<-function (W, lmin = -0.99, lmax = 0.99, incr = 0.01){
  rvec <- seq(lmin, lmax, incr)
  iN <- Matrix::Matrix(diag(nrow(W)), sparse = TRUE)
  res <- NULL
  res$rho <- NULL
  res$lndet <- NULL
  if(!is(W,"sparseMatrix")) {ww<-Matrix::Matrix(W, sparse = TRUE)} else {ww<-W}
  
  # cm <- chol(diag(nrow(W))-0.1*W, pivot=TRUE)
  # p <- attr(cm,"pivot")
  
  for (i in 1:length(rvec)) {
    rho <- rvec[i]
    z <- iN - rho * ww
    #z2<-Matrix::expand(Matrix::lu(z[,p]))
    z2<-slot(Matrix::lu(z),"U")
    res$rho[i] <- rho
    #res$lndet[i] <- sum(log(abs(diag(as.matrix(z2$U)))))
    res$lndet[i] <- sum(log(abs(diag(z2))))
  }
  return(res)
}


lndetmc <- function(W,lmin=0,lmax=0.99,m=30,p=30,sd=12345,incr=0.01){
  n<-nrow(W)
  alpha<-seq(lmin,lmax,incr)
  set.seed(sd)
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  v  <- vector(mode="list", length=m)
  trc <- sum(diag(W %*% W))
  c<-x
  for (k in 1:m) {
    c <- W %*% c
    if(k==1) {   v[[1]]<-rep(0, p)
    }else if(k==2){   v[[2]]<-rep(trc,p)}
    else{  v[[k]] <- colSums(x*c)   }
  }
  v2<-colSums(x*x)

  f1<-function(rho,m,p,v,v2,n){
    vv <- rep(0,p)
    for (k in 1:m) {   vv <- -n*(rho^k)*(v[[k]]/k) + vv   }
    return(mean(vv/v2))
  }

  lndet<-rep(0,length(alpha))
  for(i in 1:length(alpha)){ lndet[i]<-f1(alpha[i],m,p,v,v2,n) }
  result<-list()
  result$rho<-alpha
  result$lndet<-lndet
  return(result)

}
