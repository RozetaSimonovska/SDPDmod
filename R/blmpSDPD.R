#' @name blmpSDPD
#' @title Bayesian log-marginal posterior probabilities for spatial panel models
#'
#' @description Calculates log-marginal posterior probabilities for model comparison purposes.
#'
#' @param formula a symbolic description for the model to be estimated
#' @param data a data.frame
#' @param W spatial weights matrix (row-normalized)
#' @param index the indexes (names of the variables for the spatial and time component)
#' @param model a list of models for which the Bayesian log-marginal posterior probabilities need to be calculated, list("ols","slx","sar","sdm","sem","sdem")
#' @param effect type of fixed effects, c("none","individual","time","twoways"), default ="none"
#' @param ldet Type of computation of log-determinant, c("full","mc"). Default "full" for smaller problems, "mc" for large problems.
#' @param lndetspec specifications for the calculation of the log-determinant
#' @param dynamic logical, if TRUE time lag of the dependent variable is included. Default = FALSE
#' @param tlaginfo specification for the time lag, default = list(ind=NULL), \emph{ind} - i-th column in the data frame which represents the time lag
#' @param LYtrans logical, default FALSE. If Lee-Yu transformation should be used for demeaning of the variables
#' @param incr increment for vector of values for rho
#' @param rintrv logical, default TRUE, calculates eigenvalues of W. If FALSE, the interval for rho is (-1,1).
#' @param prior type of prior to be used c("uniform","beta"). Default "uniform"
#' @param bprarg argument for the beta prior. Default = 1.01
#'
#' @details
#' For the Spatial Durbin Error Model (SDEM) the marginal distribution is:
#' \deqn{p(\lambda|y) = \frac{1}{p(y)} p(\lambda) \Gamma(a) (2\pi)^{-a} \frac{|P|^{T-1}}{|Z'Z|^{1/2}} (e'e)^{-a}}
#' For the Spatial Durbin Model (SDM) the marginal distribution is:
#' \deqn{p(\rho|y) = \frac{1}{p(y)} p(\rho) \Gamma(a) (2\pi)^{-a} \frac{|P|}{|Z'Z|^{1/2}} (e'e)^{-a}}
#' where \eqn{p(\lambda)} is prior on \eqn{\lambda} and \eqn{p(\rho)} is prior on \eqn{\rho},
#' either uniform \eqn{\frac{1}{D}}, \eqn{D = 1/\omega_{max}-1/\omega_{min}} or beta prior;
#' No priors on beta and sige;
#' \eqn{\omega_{max}} and \eqn{\omega_{min}} are the maximum and minimum eigenvalues of
#' \emph{W} - spatial weights matrix;
#' \eqn{Z = X} for lag or error model and \eqn{Z = [X WX]} for Durbin model;
#' X - matrix of \eqn{k} covariates.
#'
#' More details, see LeSage (2014).
#'
#' Based on MatLab function log_marginal_panelprob.m.
#'
#' @returns  A list
#' \item{lmarginal}{log-marginal posterior}
#' \item{probs}{model probability}
#'
#' @author Rozeta Simonovska
#'
#' @import plm
#' @import RSpectra
#' @import Matrix
#'
#' @references
#' LeSage, J. P., & Parent, O. (2007). Bayesian model averaging for spatial econometric models. \emph{Geographical Analysis, 39(3)}, 241-267.
#'
#' LeSage, J. P. (2014). Spatial econometric panel data model specification: A Bayesian approach. \emph{Spatial Statistics, 9}, 122-145.
#'
#'@examples
#'\donttest{
#'## US States Production data
#'data(Produc, package = "plm")
#'## Spatial weights row-normalized matrix of 48 US states
#'data(usaww, package = "splm")
#'isrownor(usaww)
#'form1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#'res1  <- blmpSDPD(formula = form1, data=Produc, W = usaww,
#'                  index = c("state","year"),
#'                  model = list("sar","sdm","sem","sdem"),
#'                  effect = "twoways")
#'res1
#'res2  <- blmpSDPD(formula = form1, data = Produc, W = usaww,
#'                  index = c("state","year"),
#'                  model = list("sar","sdm","sem","sdem"),
#'                  effect = "twoways", dynamic = TRUE)
#'res2}
#'
#' @export

blmpSDPD<-function(formula, data, W, index, model, effect,
                   ldet = NULL, lndetspec = list(m=NULL,p=NULL,sd=NULL),
                  dynamic = FALSE, tlaginfo = list(ind = NULL),
                  LYtrans = FALSE, incr = NULL, rintrv = TRUE,
                  prior="uniform", bprarg = 1.01){
  mod_nam<-c("ols","sar","sdm","sem","sdem","slx")
  if(!any(model %in% mod_nam)) {
    stop(paste0('Wrong value for model! ',
                'Enter at least one of the following values ',
                'list("ols","sar","sdm","sem","sdem","slx")'))  }

  if(!inherits(formula, "formula")) stop("Error in formula!")
  if(is.null(index) || length(index)!=2){
    stop("Index is missing or error in index!")}

  pmod <- plm::plm(formula, data, index = index)
  ypom <- data.matrix(pmod$model[,1:2])
  dep.name <- colnames(ypom)[1]
  X <- matrix(data.matrix(pmod$model)[,-1],
              ncol = length(colnames(pmod$model))-1)
  cov.names<-colnames(pmod$model)[-1]
  y <- ypom[,1]
  sind <- attr(pmod$model, "index")[, 1]
  tind <- attr(pmod$model, "index")[, 2]
  oo <- order(tind, sind)
  X <- X[oo, , drop = FALSE]
  y <- y[oo]
  sind <- sind[oo]
  tind <- tind[oo]
  n <- length(unique(sind))  ##number of individuals
  k <- dim(X)[[2]]  ##number of covariates
  t <- max(tapply(X[, 1], sind, length))  ##number of years

  if (nrow(W) != n) stop("Non conformable spatial weights")

  if (!is.matrix(W)) {
    if ("listw" %in% class(W)) { W <- listw2mat(W)
    }else { stop("W has to be either a 'matrix' or a 'listw' object!")  }
  }

  balanced <- plm::pdim(pmod)$balanced
  if (!balanced) stop("Estimation method unavailable for unbalanced panels!")

  if(is.null(effect)){ effect<-"none"
  }else if(!is.null(effect) & !(effect %in%
                                c("none","individual","time","twoways"))){
    stop("Wrong fixed effect entered!")}

  if(dynamic){
    if(!is.null(tlaginfo$ind)){
      if(is.numeric(tlaginfo$ind)){
        tlgy<-data[,tlaginfo$ind]
        tlgy<-tlgy[oo]
        for(i in 1:(t-1)){
          for(j in 1:n){
            if(tlgy[n*(i)+j]!=y[n*(i-1)+j]) stop("Wrong index for time lag!")
          }
        }
      } else{ stop("Non numeric index for time lag of the dependent variable!")}
    }else{
      tlgy<-y[1:(n*(t-1))]
      X<-X[(n+1):(n*t),]
      rownames(X)<-seq(1,nrow(X),1)
      y<-y[(n+1):(n*t)]
      y<-t(t(y))
      rownames(y)<-seq(1,nrow(y),1)
      sind1<-sind; tind1<-tind; tind<-rep(NA,n*(t-1))
      sind<-sind1[seq(n+1,n*t,1)]; levels(sind)<-as.character(sind)
      tind<-tind1[seq(n+1,n*t,1)]; levels(tind)[1]<-NA
      t<-t-1
    }
    X<-cbind(tlgy,X)  ###add time lag to X matrix
    k<-k+1
  }

  wrnor<-isrownor(W)

  ####Demeaning method
  if(effect=="none"){
    message("\nNo demeaning used.\n")
  }else if(effect %in% c("individual","time")){
    if(LYtrans & wrnor & !dynamic){
      re2<-demeanF(y,X,n,t,effect,W)
      y<-re2$yf; X<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
    }else{
      re1<-demean(y,X,n,t,effect,sind,tind)
      y<-re1$yw; X<-re1$xw
    }
  }else if(effect %in% c("twoways")){
    if(LYtrans & dynamic & wrnor){
      sind2<-sind[-seq(1,length(sind),n)]; levels(sind2)[1]<-NA
      tind2<-tind[-seq(1,length(tind),n)]; levels(tind2)<-as.character(tind2)
      re2<-demeanF(y,x=X,n,t,effect="time",W)
      yy<-re2$yf; Xx<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
      re1<-demean(y=yy,x=Xx,n,t,effect="individual",sind2,tind2)
      y<-re1$yw; X<-re1$xw
    } else if(LYtrans & !dynamic & wrnor){
      re2<-demeanF(y,X,n,t,effect,W)
      y<-re2$yf; X<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
    }else{
      LYtrans<-FALSE
      re1<-demean(y,X,n,t,effect,sind,tind)
      y<-re1$yw; X<-re1$xw
    }
  }

  ####Adjustment for degrees of freedom due to included fixed effects
  if(effect=="none"){
    dofadj<-1   ## intercept will be included
  }else if(effect=="individual"){
    dofadj<-n   ## correction for n spatial fixed effects
  }else if(effect=="time"){
    dofadj<-t  ## correction for t time-period fixed effects
  }else if(effect=="twoways"){
    dofadj<-n+t-1   ## correction for spatial and time-period fixed effects
  }else stop("wrong entry in  fixed effects")


  ####increment
  if(is.null(incr) & n<500){incr <- 0.001
  }else if(is.null(incr) & n>=500){ incr <- 0.01 }

  ##eigenvalues
  if(rintrv & prior=="uniform"){
    ei.max <- Re(RSpectra::eigs(W,1,which = "LR")$values)
    ei.min <- Re(RSpectra::eigs(W,1,which = "SR")$values)
    if(length(ei.min)==0){
      warning("Minimun eigen value not found."); ei.min<-(-1)}
    rmin <- 1/ei.min + incr;    rmax <- 1/ei.max - incr
  } else { rmin <- (-1) + incr;     rmax <- 1 - incr }

  #####Log-determinant calculation
  if(is.null(ldet)){
    if(n<1000){
      out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
    } else {
      if(!is.null(lndetspec$p) & !is.null(lndetspec$m) &
         !is.null(lndetspec$sd)) {
        rmin<-0
        out <- lndetmc(W,lmin=rmin,lmax=rmax,p=lndetspec$p,
                       m=lndetspec$m,sd=lndetspec$sd,incr)
      }else {
        rmin<-0
        out <- lndetmc(W,lmin=rmin,lmax=rmax,m=30,p=30,sd=12345,incr)
      }
    }
  } else if(ldet=="full"){
    out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
  } else if(ldet=="mc"){
    if(!is.null(lndetspec$p) & !is.null(lndetspec$m) &
       !is.null(lndetspec$sd)) {
      rmin<-0
      out <- lndetmc(W,lmin=rmin,lmax=rmax,p=lndetspec$p,
                     m=lndetspec$m,sd=lndetspec$sd,incr)
    }else {
      rmin<-0
      out <- lndetmc(W,lmin=rmin,lmax=rmax,m=30,p=30,sd=12345,incr)
    }
  } else{
    out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
    warning(paste0("Wrong entry for log-determinant. ",
                   "Continuing with calculation of lndetfull!"))
  }

  ####interpolation
  if(incr>0.001){
    rvect <- seq(rmin,rmax,0.001)
    outi<-spline(x = out$rho, y = out$lndet, n = length(rvect),
                 xmin = min(rvect), xmax = max(rvect), method = "fmm")
    detval <-cbind(outi$x,outi$y)
  }else{  detval <-cbind(out$rho,out$lndet)}

  ######################
  Wnt<-kronecker(diag(t),W)
  Wx<-as.matrix(Wnt%*%X)
  Wy<-as.matrix(Wnt%*%y)
  WWx<-as.matrix(Wnt%*%Wx)

  #############################
  lmarg_df <- as.data.frame(matrix(NA, ncol=length(mod_nam)))
  colnames(lmarg_df) <- mod_nam

  ###############
  if(effect=="none") {
    xxm <- cbind(rep(1,n*t),X); xxdm <- cbind(rep(1,n*t),X,Wx)
  }else {  xxm<- X; xxdm<- cbind(X,Wx) }

  ####Prior on rho
  if(prior=="uniform"){
    D <- (1/ei.max - 1/ei.min)
    bpr <- 1
  }else if(prior=="beta"){
    bpr <- betapr(detval[,1],bprarg)
    D <- 1
  }else { stop("Error in type of prior!")}

  #################################################################
  ######log-marginal
  if(any(model %in% "ols")){ # log-marginal for least-squares case
    xx <- xxm
    dof <- (n*t-k-dofadj)/2
    xpx <- t(xx)%*%xx
    lndetx <- log(det(xpx))
    logC <- lgamma(dof) - dof*log(2*pi) - lndetx/2
    bhat <- solve(xpx)%*%(t(xx)%*%y)
    epe <- c(t(y-xx%*%bhat)%*%(y-xx%*%bhat))
    logm_out <- c(logC - dof*log(epe))
    lmarg_df$ols<-logm_out
  }

  if(any(model %in% "slx")){
    xx <- xxdm
    dof <- (n*t-2*k-dofadj)/2
    xpx <- t(xx)%*%xx
    lndetx <- log(det(xpx))
    logC <- lgamma(dof) - dof*log(2*pi) - lndetx/2
    bhat <- solve(xpx)%*%(t(xx)%*%y)
    epe <- c(t(y-xx%*%bhat)%*%(y-xx%*%bhat))
    logm_out <- c(logC - dof*log(epe))
    lmarg_df$slx <- logm_out
  }

  if(any(model %in% "sar")){
    xx<-xxm
    dof <- (n*t-k-dofadj)/2
    xpx <- t(xx)%*%xx
    lndetx <- log(det(xpx))
    logC <- (-log(D)) + lgamma(dof) - dof*log(2*pi) - lndetx/2
    ###############
    bo <- solve(xpx)%*%(t(xx)%*%y)
    bd <- solve(xpx)%*%(t(xx)%*%Wy)
    eo <- y - xx%*%bo
    ed <- Wy - xx%*%bd
    epeo <- as.vector(t(eo)%*%eo)
    eped <- as.vector(t(ed)%*%ed)
    epeod <- as.vector(t(ed)%*%eo)
    Q1 <- epeo - 2*detval[,1]*epeod + (detval[,1]*detval[,1])*eped
    Q2 <- detval[,2]
    logm <-rep(0,length(detval[,1]))
    logm <- (-dof*log(Q1)) + t*Q2 + log(bpr)
    adj <- max(logm)
    madj <- logm - adj
    fint <- exp(madj)
    isum<-sum((detval[2:length(detval[,1]),1] +
                 detval[1:(length(detval[,1]) - 1),1]) * (
                   fint[2:length(detval[,1])] -
                     fint[1:(length(detval[,1]) - 1)])/2)
    ############
    logm_out <- isum + adj + logC
    lmarg_df$sar<-logm_out
  }

  if(any(model %in% "sdm")){
    xx<-xxdm
    dof <- (n*t-2*k-dofadj)/2
    xpx <- t(xx)%*%xx
    lndetx <- log(det(xpx))
    logC <- (-log(D)) + lgamma(dof) - dof*log(2*pi) - lndetx/2
    #######
    bo <- solve(xpx)%*%(t(xx)%*%y)
    bd <- solve(xpx)%*%(t(xx)%*%Wy)
    eo <- y - xx%*%bo
    ed <- Wy - xx%*%bd
    epeo <- as.vector(t(eo)%*%eo)
    eped <- as.vector(t(ed)%*%ed)
    epeod <- as.vector(t(ed)%*%eo)
    Q1 <- epeo - 2*detval[,1]*epeod + (detval[,1]*detval[,1])*eped
    Q2 <- detval[,2]
    logm <- rep(0,length(detval[,1]))
    logm <- (-dof*log(Q1)) + t*Q2 + log(bpr)
    adj <- max(logm)
    madj <- logm - adj
    fint <- exp(madj)
    isum <- sum((detval[2:length(detval[,1]),1] +
                   detval[1:(length(detval[,1]) - 1),1]) * (
                     fint[2:length(detval[,1])] -
                       fint[1:(length(detval[,1]) - 1)])/2)
    #######
    logm_out <- isum + adj + logC
    lmarg_df$sdm <- logm_out
  }
  Wx0<-Wx
  if(any(model %in% "sem")){
    xx <- xxm
    dof <- (n*t-k-dofadj)/2
    logC <- (-log(D)) + lgamma(dof) - dof*log(2*pi)

    if(effect=="none"){Wx<-cbind(rep(1,n*t),Wx)}
    ######
    xpx <- t(xx)%*%xx
    xpWx <- t(xx)%*%Wx
    xpWpx <- t(Wx)%*%xx
    xpWpWx <- t(Wx)%*%Wx
    xpy <- t(xx)%*%y
    xpWy <- t(xx)%*%Wy
    xpWpy <- t(Wx)%*%y
    xpWpWy <- t(Wx)%*%Wy
    ypy <- as.numeric(crossprod(as.vector(y)))
    ypWy <- as.numeric(t(as.vector(y))%*%Wy)
    ypWpy <- as.numeric(t(Wy)%*%as.vector(y))
    ypWpWy <- as.numeric(t(Wy)%*%Wy)
    Q1 <- rep(0,length(detval[,1]))
    Q3 <- rep(0,length(detval[,1]))
    for(it in 1:length(detval[,1])){
      rho <- detval[it,1]
      Axx <- xpx - rho*xpWx - rho*xpWpx + rho*rho*xpWpWx
      Q3[it] <- log(det(Axx))
      Axy <- xpy - rho*xpWy - rho*xpWpy + rho*rho*xpWpWy
      Ayy <- ypy - rho*ypWy - rho*ypWpy + rho*rho*ypWpWy
      b <- solve(Axx)%*%Axy
      Q1[it] <- as.numeric((Ayy - t(b)%*%Axx%*%b))
    }
    Q2 <- detval[,2]
    logm <-rep(0,length(detval[,1]))
    logm <- - dof*log(Q1) + t*Q2 - Q3/2 + log(bpr)
    adj <- max(logm)
    madj <- logm - adj
    fint <- exp(madj)
    isum <- sum((detval[2:length(detval[,1]),1] +
                   detval[1:(length(detval[,1]) - 1),1]) * (
                     fint[2:length(detval[,1])] -
                       fint[1:(length(detval[,1]) - 1)])/2)
    ###########
    logm_out <- isum + adj + logC
    lmarg_df$sem <- logm_out
  }
  Wx<-Wx0
  if(any(model %in% "sdem")){
    xx <- xxdm
    if(effect=="none"){Wxsem<-cbind(rep(1,n*t),Wx,WWx)
    }else { Wxsem <- cbind(Wx,WWx)  }
    dof <- (n*t-2*k-dofadj)/2
    logC <- (-log(D)) + lgamma(dof) - dof*log(2*pi)
    #####
    xpx <- t(xx)%*%xx
    xpWx <- t(xx)%*%Wxsem
    xpWpx <- t(Wxsem)%*%xx
    xpWpWx <- t(Wxsem)%*%Wxsem
    xpy <- t(xx)%*%y
    xpWy <- t(xx)%*%Wy
    xpWpy <- t(Wxsem)%*%y
    xpWpWy <- t(Wxsem)%*%Wy
    ypy <- as.numeric(crossprod(as.vector(y)))
    ypWy <- as.numeric(t(as.vector(y))%*%Wy)
    ypWpy <- as.numeric(t(Wy)%*%as.vector(y))
    ypWpWy <- as.numeric(t(Wy)%*%Wy)
    Q1 <- rep(0,length(detval[,1]))
    Q3 <- rep(0,length(detval[,1]))
    for(it in 1:length(detval[,1])){
      rho <- detval[it,1]
      Axx <- xpx - rho*xpWx - rho*xpWpx + rho*rho*xpWpWx
      Q3[it] <- log(det(Axx))
      Axy <- xpy - rho*xpWy - rho*xpWpy + rho*rho*xpWpWy
      Ayy <- ypy - rho*ypWy - rho*ypWpy + rho*rho*ypWpWy
      b <- solve(Axx)%*%Axy
      Q1[it] <- as.numeric((Ayy - t(b)%*%Axx%*%b))
    }
    Q2 <- detval[,2]
    logm <-rep(0,length(detval[,1]))
    logm <- (- dof*log(Q1)) + t*Q2 - Q3/2 + log(bpr)
    adj <- max(logm)
    madj <- logm - adj
    fint <- exp(madj)
    isum<-sum((detval[2:length(detval[,1]),1] +
                 detval[1:(length(detval[,1]) - 1),1]) * (
                   fint[2:length(detval[,1])] -
                     fint[1:(length(detval[,1]) - 1)])/2)
    ####
    logm_out <- isum + adj + logC
    lmarg_df$sdem<-logm_out
  }

  lmarginal<-lmarg_df[which(!is.na(lmarg_df))]
  adj <- max(lmarginal)
  madj <- lmarginal - adj
  xxv <- exp(madj)
  psum <- sum(xxv)
  probs <- xxv/psum

  result<-list(lmarginal,probs)
  names(result)<-c("lmarginal","probs")
  return(result)
}
