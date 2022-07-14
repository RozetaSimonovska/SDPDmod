#####fixed effects function
feffects<-function(rho,beta,sige,W0,y,X,n0,t0,y0,X0,mny,mnx,mty,mtx,effect,tind,sind,Wy0)
{
  mnWy<-rep(0,n0)
  mnWy<-tapply(Wy0,sind,mean)
  mtWy<-rep(0,t0)
  mtWy<-tapply(Wy0,tind,mean)

  intercept<-mean(y0)-mean(Wy0)*rho-colMeans(X0)%*%beta
  int.se<-sqrt(sige/(n0*t0)+sige*t(colMeans(X0))%*%solve(t(X)%*%X)%*%(colMeans(X0))) ## s.e. itercept
  int.t<-intercept/int.se ## t stat intercept
  int.p<- 2*pnorm(abs(int.t),lower.tail=FALSE)
  int.tab <- cbind(intercept,int.se,int.t,int.p)
  colnames(int.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  rownames(int.tab) <- "(Intercept)"
  sfe<-NA; tfe<-NA; tfe.se<-NA; sfe.se<-NA

  result<-list()
  result$int.tab<-int.tab

  if(effect=="individual"){
    sfe<-mny-mnWy*rho-c(mnx%*%beta)-rep(intercept,n0) ## spatial fixed effects coefficients
    xhat<-X0%*%beta+rep(sfe,t0)+rep(intercept,n0*t0)
    sfe.se<-sqrt(sige/t0*rep(1,n0)+diag(sige*mnx%*%solve(t(X)%*%X)%*%t(mnx))) ###standard errors
    sfe.t<-sfe/sfe.se  ###t statistic
    sfe.p<- 2*pnorm(abs(sfe.t),lower.tail=FALSE)
    sfe.tab <- cbind(sfe,sfe.se,sfe.t,sfe.p)
    colnames(sfe.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    result$sfe.tab<-sfe.tab
  } else if(effect=="time"){
    tfe<-mty-mtWy*rho-c(mtx%*%beta)-rep(intercept,t0)
    xhat<-X0%*%beta+rep(tfe,each=n0)+rep(intercept,n0*t0)
    tfe.se<-sqrt(sige/n0*rep(1,t0)+diag(sige*mtx%*%solve(t(X)%*%X)%*%t(mtx)))
    tfe.t<-tfe/tfe.se
    tfe.p<- 2*pnorm(abs(tfe.t),lower.tail=FALSE)
    tfe.tab <- cbind(tfe,tfe.se,tfe.t,tfe.p)
    colnames(tfe.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    result$tfe.tab<-tfe.tab
  } else if(effect=="twoways"){
    sfe<-mny-mnWy*rho-c(mnx%*%beta)-rep(intercept,n0)
    tfe<-mty-mtWy*rho-c(mtx%*%beta)-rep(intercept,t0)
    sfe.se<-sqrt(sige/t0*rep(1,n0)+diag(sige*mnx%*%solve(t(X)%*%X)%*%t(mnx))) ###standard errors
    sfe.t<-sfe/sfe.se  ###t statistic
    sfe.p<- 2*pnorm(abs(sfe.t),lower.tail=FALSE)
    sfe.tab <- cbind(sfe,sfe.se,sfe.t,sfe.p)
    colnames(sfe.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    result$sfe.tab<-sfe.tab

    tfe.se<-sqrt(sige/n0*rep(1,t0)+diag(sige*mtx%*%solve(t(X)%*%X)%*%t(mtx)))
    tfe.t<-tfe/tfe.se
    tfe.p<- 2*pnorm(abs(tfe.t),lower.tail=FALSE)
    tfe.tab <- cbind(tfe,tfe.se,tfe.t,tfe.p)
    colnames(tfe.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    result$tfe.tab<-tfe.tab

    xhat<-X0%*%beta+rep(sfe,t0)+rep(tfe,each=n0)+rep(intercept,n0*t0)
  } else if(effect=="none"){
    xhat<-X0%*%beta
  }
  res.e <- y0 - xhat - rho* Wy0

  ywhat <-  X %*% beta
  res1 <- as.matrix(y - mean(y))
  res2 <- as.matrix(ywhat - mean(ywhat))
  rsq1 <- t(res1)%*%res2
  rsq2 <- t(res1)%*%res1
  rsq3 <- t(res2)%*%res2
  corr2 <- c(rsq1^2/(rsq2*rsq3))

  result$intercept<-intercept
  result$int.se<-int.se
  result$int.t<-int.t
  result$sfe<-sfe
  result$sfe.se<-sfe.se
  result$tfe<-tfe
  result$tfe.se<-tfe.se
  result$xhat<-xhat
  result$res.e<-res.e
  result$corr2<-corr2
  return(result)
}


