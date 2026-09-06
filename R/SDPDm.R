#' @name SDPDm
#' @title Spatial dynamic panel data lag model with fixed effects maximum likelihood estimation.
#'
#' @description This function estimates spatial panel model with fixed effects for static or dynamic model. It includes the transformation approach suggested by Yu et al (2008) and Lee and Yu (2010).
#'
#' @param formula a symbolic description for the (static) model to be estimated, not including the dynamic component
#' @param data a data.frame
#' @param W spatial weights matrix
#' @param index the indexes (Names of the variables for the spatial and time component. The spatial is first and the time second.)
#' @param model a models to be calculated, c("sar","sdm"), default = "sar"
#' @param effect type of fixed effects, c("none","individual","time","twoways"), default ="individual"
#' @param ldet type of computation of log-determinant, c("full","mc"). Default "full" for smaller problems, "mc" for large problems.
#' @param lndetspec specifications for the calculation of the log-determinant for mcmc calculation. Default list(p=NULL,m=NULL,sd=NULL), if the number of spatial units is >1000 then list(p=30,m=30,sd=12345)
#' @param dynamic logical, if TRUE time lag of the dependent variable is included. Default = FALSE
#' @param tlaginfo specification for the time lag, default = list(ind=NULL,tl=FALSE,stl=FALSE), see details
#' @param LYtrans logical, default FALSE. If the Lee-Yu transformation should be used for bias correction
#' @param incr increment for vector of values for rho
#' @param rintrv logical, default TRUE, calculates eigenvalues of W. If FALSE, the interval for rho is (-1,1)
#' @param demn logical, if Lee-Yu transformation for demeaning of the variables to remove fixed effects is performed (only used in static models). Default FALSE
#' @param DIRtrans logical, if direct transformation of variables should be used. Default, FALSE (only used in dynamic models with "twoways" effects)
#'
#' @details
#' Based on MatLab functions sar_jihai.m, sar_jihai_time.m and sar_panel_FE.m
#'
#'  In \emph{tlaginfo = list(ind = NULL, tl = TRUE, stl = TRUE)}:
#'
#' \emph{ind} i-th column in \emph{data} which represents the time lag, 
#' if not specified then the lag from the dependent variable is created and the 
#' panel is reduced from n*t to n*(t-1)
#'
#' \emph{tl} logical, default TRUE. If TRUE \eqn{y_{t-1}} 
#' (the lagged dependent variable in time is included)
#'
#' \emph{stl} logical, default TRUE. If TRUE \eqn{Wy_{t-1}} 
#' (the lagged dependent variable in space and time is included)
#'
#' @returns An object of class "SDPDm"
#' \item{coefficients}{coefficients estimate of the model parameters (\emph{coefficients1} for dynamic model)}
#' \item{rho}{spatial coefficient}
#' \item{sige}{residuals variance}
#' \item{llik}{the value of the log likelihood function}
#' \item{...}{}
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{vignette("spatial_model", package = "SDPDmod")}
#'
#' @import plm
#' @import RSpectra
#' @import Matrix
#' @importFrom stats optimize pchisq pnorm printCoefmat rnorm spline
#'
#' @references
#' Yu, J., De Jong, R., & Lee, L. F. (2008). Quasi-maximum likelihood estimators for spatial dynamic panel data with fixed effects when both n and T are large. \emph{Journal of Econometrics}, 146(1), 118-134.
#'
#' Lee, L. F., & Yu, J. (2010). Estimation of spatial autoregressive panel data models with fixed effects. \emph{Journal of Econometrics}, 154(2), 165-185.
#'
#' Lee, L. F., & Yu, J. (2010). A spatial dynamic panel data model with both time and individual fixed effects. \emph{Econometric Theory}, 564-597.
#'
#' @examples
#' \donttest{
#' library("SDPDmod")
#' data(Produc, package = "plm")
#' data(usaww, package = "splm")
#' form1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#' mod1  <- SDPDm(formula = form1, data = Produc, W = usaww, index = c("state","year"),
#'                model = "sar", effect = "individual", LYtrans = TRUE)
#' summary(mod1)
#' imp1  <- impactsSDPDm(mod1)
#' summary(imp1)
#' mod2  <- SDPDm(formula = form1, data = Produc, W = usaww, index = c("state","year"),
#'                model = "sdm", effect = "twoways", LYtrans = TRUE,
#'                dynamic = TRUE, tlaginfo=list(ind = NULL, tl = TRUE, stl = TRUE))
#' summary(mod2)}
#'
#' @export


SDPDm<-function(formula, data, W, index, 
                model = "sar", effect = "individual",
                ldet = NULL, lndetspec=list(p=NULL,m=NULL,sd=NULL),
                dynamic = FALSE,
                tlaginfo = list(ind = NULL,tl = TRUE,stl = TRUE),
                LYtrans = FALSE,
                incr = NULL, rintrv = TRUE,
                demn = FALSE, DIRtrans = FALSE)
{
  cl <- match.call()
  if(!inherits(formula, "formula")){
    stop("Error in formula!")}  ###static model formula
  if(is.null(index) || length(index)!=2){
    stop("Index is missing or error in index!")}

  data2 <- data[order(data[,index[1]],data[,index[2]]),]
  pmod <- plm::plm(formula, data2, index = index)
  ypom<-data.matrix(pmod$model[,1:2])
  dep.name<-colnames(ypom)[1]
  X <- matrix(data.matrix(pmod$model)[,-1],
              ncol = length(colnames(pmod$model))-1)
  cov.names<-colnames(pmod$model)[-1]
  y <- ypom[,1]
  att_ind <- attr(pmod$model, "index")
  sind <- att_ind[, 1]  ##index individuals/regions
  tind <- att_ind[, 2]  ##index time
  oo <- order(tind, sind)   ####order obs for each region in a year
  X <- X[oo, , drop = FALSE]

  y <- y[oo]
  sind <- sind[oo]
  tind <- tind[oo]
  n <- length(unique(sind))  ##number of individuals
  k <- dim(X)[[2]]  ##number of covariates
  t <- max(tapply(X[, 1], sind, length))  ##number of years

  ###check if number of rows of wight matrix is equal to number of
  ###individuals/regions
  if (nrow(W) != n) stop("Non conformable spatial weights!")

  if (!is.matrix(W)) {
    if ("listw" %in% class(W)) { W <- listw2mat(W)
    }else { stop("W has to be either a 'matrix' or a 'listw' object!")
    }
  }

  balanced <- plm::pdim(pmod)$balanced
  ###stop if unbalanced panel
  if (!balanced) stop("Estimation method unavailable for unbalanced panels!")

  if(is.null(effect)){ effect<-"individual"
  }else if(!is.null(effect) & !(effect %in%
                                c("none","individual","time","twoways"))) {
    stop("Wrong fixed effects entered!")}
  
  if(!is.null(model)){
    if(!(model %in% c("sar","sdm"))){
      stop("Wrong model entered! 
           Enter 'sar' for spatial autoregressive model or 'sdm' for spatial Durbin model!")
    }
  }

  if(dynamic){
    if(!is.null(tlaginfo$ind)){
      if(is.numeric(tlaginfo$ind)){
        tlagy0 <- matrix(data2[oo,tlaginfo$ind])
        tlagy <- tlagy0[,1]
        for(i in 1:(t-1)){
          for(j in 1:n){
            if(tlagy[n*(i)+j]!=y[n*(i-1)+j])   stop("Wrong index for time lag!")
          }
        }
      } else {stop("Non numeric index for time lag of the dependent variable!")}
    }else{
      tlagy<-y[1:(n*(t-1))]
      X<-as.matrix(X[(n+1):(n*t),],ncol=k)
      rownames(X)<-seq(1,nrow(X),1)
      y<-y[(n+1):(n*t)]
      y<-t(t(y))
      rownames(y)<-seq(1,nrow(y),1)
      sind1<-sind; tind1<-tind; tind<-rep(NA,n*(t-1))
      sind<-sind1[seq(n+1,n*t,1)]; levels(sind)<-as.character(sind)
      tind<-tind1[seq(n+1,n*t,1)]; levels(tind)[1]<-NA
      t<-t-1
    }
    X<-cbind(X,tlagy)  ###add time lag to X matrix
    k<-k+1
  }

  wrnor<-ifelse(isrownor(W),TRUE,FALSE)

  Wnt<-kronecker(diag(t),W)

  if(dynamic){
    tlagy <- X[,ncol(X)]; X <- X[,-ncol(X)]; k <- k-1
    if(tlaginfo$stl){ Wty <- as.matrix(Wnt%*%tlagy)   }
    Wx <- as.matrix(Wnt%*%X);  Wy <- as.matrix(Wnt%*%y)
    if(tlaginfo$tl & tlaginfo$stl){
      if(is.null(model) || model=="sar"){     Z<-cbind(tlagy,Wty,X)
      }else if(model=="sdm"){                 Z<-cbind(tlagy,Wty,X,Wx)
      } else {stop("Wrong model!")}
    } else if(!tlaginfo$tl & tlaginfo$stl){
      if(is.null(model) || model=="sar"){     Z<-cbind(Wty,X)
      }else if(model=="sdm"){                 Z<-cbind(Wty,X,Wx)
      } else {stop("Wrong model!")}
    } else if(tlaginfo$tl & !tlaginfo$stl){
      if(is.null(model) || model=="sar"){     Z<-cbind(tlagy,X)
      }else if(model=="sdm"){                 Z<-cbind(tlagy,X,Wx)
      } else {stop("Wrong model!")}
    } else { stop("Wrong values for dynamic spatial and time interactions!")
    }
  } else {
    Wx <- as.matrix(Wnt%*%X);    Wy <- as.matrix(Wnt%*%y)
    if(is.null(model) || model=="sar"){       Z<-X
    }else if(model=="sdm"){                   Z<-cbind(X,Wx)
    } else {stop("Wrong model")}
  }

  y0<-y; Z0<-Z; n0<-n; t0<-t; W0<-W; Wy0<-Wy; sind0<-sind; tind0<-tind
  kz<-ncol(Z)

  if(dynamic & LYtrans & DIRtrans & effect=="twoways"){
    stop("LYtrans and DIRtrans can not be used at the same time!")}

  ####Demeaning
  if(effect=="none"){
    message("No demeaning used.")
  }else if(effect %in% c("individual","time")){
    if(LYtrans & wrnor & demn & !dynamic){
      re2<-demeanF(y,x=Z,n,t,effect,W)
      y<-re2$yf; x<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
    }else{
      re1<-demean(y,x=Z,n,t,effect,sind,tind)
      y<-re1$yw; x<-re1$xw; mny<-re1$mny; mnx<-re1$mnx
      mty<-re1$mty; mtx<-re1$mtx
      demn<-FALSE
    }
  }else if(effect %in% c("twoways")){
    if(LYtrans & dynamic & wrnor){
      sind2<-sind[-seq(1,length(sind),n)]; levels(sind2)[1]<-NA
      tind2<-tind[-seq(1,length(tind),n)]; levels(tind2)<-as.character(tind2)
      re2<-demeanF(y,x=Z,n,t,effect="time",W)
      yy<-re2$yf; Xx<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
      re1<-demean(y=yy,x=Xx,n,t,effect="individual",sind2,tind2)
      y<-re1$yw; x<-re1$xw;  mny<-re1$mny; mnx<-re1$mnx
    }else if(!LYtrans & DIRtrans & dynamic){
      Q<-kronecker(diag(t)-matrix(1/t, nrow = t, ncol = t),diag(n)-
                     matrix(1/n, nrow = n, ncol = n))
      y<-Q%*%y
      x<-Q%*%Z
    }else if(LYtrans & !dynamic & wrnor & demn){
      re2<-demeanF(y,x=Z,n,t,effect,W)
      y<-re2$yf; x<-re2$xf; W<-re2$Wf; n<-re2$nv; t<-re2$tv
    }else{
      demn<-FALSE
      re1<-demean(y,x=Z,n,t,effect,sind,tind)
      y<-re1$yw; x<-re1$xw; mny<-re1$mny; mnx<-re1$mnx
      mty<-re1$mty; mtx<-re1$mtx
    }
  }

  if(effect!="none") Z<-x
  Wnt<-kronecker(diag(t),W);  Wy <- as.matrix(Wnt%*%y)

  ####increment
  if(is.null(incr) & n<500){incr <- 0.001
  } else if(is.null(incr) & n>=500){ incr <- 0.01 }

  ####eigenvalue calculation
  if(rintrv){
    ei.max <- Re(RSpectra::eigs(W,1,which = "LR")$values)
    ei.min <- Re(RSpectra::eigs(W,1,which = "SR")$values)
    if(length(ei.min)==0){
      warning("Minimun eigen value not found."); ei.min<-(-1)}
    rmin <- 1/ei.min + incr;    rmax <- 1/ei.max - incr
  } else if(dynamic & LYtrans & effect=="twoways"){
    rmin <- 0 + incr;     rmax <- 1 - incr
  }else { rmin <- (-1) + incr;     rmax <- 1 - incr }

  ###logdet
  if(is.null(ldet)){
    if(n<1000){
      out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
    } else {
      if(!is.null(lndetspec$p) & !is.null(lndetspec$m) &
         !is.null(lndetspec$sd)) {
        out <- lndetmc(W,lmin=rmin,lmax=rmax,
                       p=lndetspec$p,m=lndetspec$m,
                       sd=lndetspec$sd,incr)
      }else {
        out <- lndetmc(W,lmin=rmin,lmax=rmax,m=30,p=30,sd=12345,incr)
      }
    }
  } else if(ldet=="full"){
    out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
  } else if(ldet=="mc"){
    if(!is.null(lndetspec$p) & !is.null(lndetspec$m) & !is.null(lndetspec$sd)) {
      out <- lndetmc(W,lmin=rmin,lmax=rmax,
                     p=lndetspec$p,m=lndetspec$m,sd=lndetspec$sd,incr)
    }else {
      out <- lndetmc(W,lmin=rmin,lmax=rmax,m=30,p=30,sd=12345,incr)
    }
  } else{
    out <- lndetfull(W,lmin=rmin,lmax=rmax,incr)
    warning(paste0("Wrong entry for log-determinant. ",
                   "Continuing with calculation of lndetfull!"))
  }

  if(incr>0.001){
    rvect <- seq(rmin,rmax,0.001)
    outi<-spline(x=out$rho, y=out$lndet, n=length(rvect),
                 xmin = min(rvect),xmax = max(rvect), method = "fmm")
    detval <-cbind(outi$x,outi$y)
  }else{  detval <-cbind(out$rho,out$lndet)}

  AI <- t(Z)%*%Z
  bo <- solve(AI)%*%(t(Z)%*%y)
  bd <- solve(AI)%*%(t(Z)%*%Wy)
  eo <- y - Z%*%bo
  ed <- as.matrix(Wy - Z%*%bd)
  epeo <- as.vector(t(eo)%*%eo)
  eped <- as.vector(t(ed)%*%ed)
  epeod <- as.vector(t(ed)%*%eo)

  optres <- optimize(f_sar,lower=detval[1,1],upper=detval[nrow(detval),1],
                     maximum = FALSE,detval=detval,
                     epeo=epeo,eped=eped,epeod=epeod,
                     n=n,t=t,dynamic=dynamic)
  rho <- optres$minimum;     liktmp <- optres$objective

  #####
  bhat <-bo - rho*bd
  res.e <- (eo - rho*ed)
  fit <- y - res.e
  yhat <- Matrix::solve(Matrix::Matrix(diag(n*t) -
                                       rho*Wnt,sparse = TRUE))%*%(Z%*%bhat)
  yhat<-c(as.matrix(yhat))
  resid <- y-yhat

  sige <-as.vector(((t(res.e)%*%(res.e))/(n*t)))
  if(LYtrans & !demn & effect == "time" & !dynamic){
    sige <- (n/(n-1)) * as.numeric(sige)  }
  if(LYtrans & !demn & effect == "individual" & !dynamic){
    sige <- (t/(t-1)) * as.numeric(sige) }

  names(sige)<-"sige";  names(rho)<-"rho"
  if(dynamic){
    residr<-as.vector(y-rho*Wy-Z%*%bhat)

    if(tlaginfo$tl & tlaginfo$stl){
      names(bhat)[1]<-paste0(dep.name,"(t-1)")
      names(bhat)[2]<-paste0("W*",dep.name,"(t-1)")
      names(bhat)[3:(k+2)]<-cov.names
      if(model=="sdm"){
        rownames(bhat)[(k+3):length(bhat)]<-paste0("W*",cov.names)
        names(bhat)[(k+3):length(bhat)]<-paste0("W*",cov.names)  }
    } else if(!tlaginfo$tl & tlaginfo$stl){
      names(bhat)[1]<-paste0("W*",dep.name,"(t-1)")
      names(bhat)[2:(k+1)]<-cov.names
      if(model=="sdm"){ row
        rownames(bhat)[(k+2):length(bhat)]<-paste0("W*",cov.names)
        names(bhat)[(k+2):length(bhat)]<-paste0("W*",cov.names)  }
    } else if(tlaginfo$tl & !tlaginfo$stl){
      names(bhat)[1]<-paste0(dep.name,"(t-1)")
      names(bhat)[2:(k+1)]<-cov.names
      if(model=="sdm"){
        rownames(bhat)[(k+2):length(bhat)]<-paste0("W*",cov.names)
        names(bhat)[(k+2):length(bhat)]<-paste0("W*",cov.names)  }
    }
  }else{
    names(bhat)[1:k]<-cov.names
    if(model=="sdm"){
      rownames(bhat)[(k+1):length(bhat)]<-paste0("W*",cov.names)
      names(bhat)[(k+1):length(bhat)]<-paste0("W*",cov.names)  }
  }

  if(LYtrans & !demn & effect == "time" & !dynamic){
    likl <- f2_sar(rho,bhat,y,Z,Wy,detval,n-1,t,sige)
  }else if(LYtrans & !demn & effect == "individual" & !dynamic){
    likl <- f2_sar(rho,bhat,y,Z,Wy,detval,n,t-1,sige)
  }else if((dynamic & LYtrans & effect %in% c("individual","twoways")) ||
           (dynamic & DIRtrans & effect %in% c("twoways"))){
    likl <-f2_sar_dyn(rho,bhat,y0,Z0,detval,n0,t0,sige,sind0,tind0,Wy0)
  }else { likl <-f2_sar(rho,bhat,y,Z,Wy,detval,n,t,sige)   }

  ###
  fsig<-f_SIG(rho,bhat,sige,W,Z,n,t,kz)
  SIG<-fsig$SIG; Gn<-fsig$Gn; Sni<-fsig$Sni
  if((LYtrans & !demn & effect=="twoways") || dynamic) {SIG <- SIG/(n*t)}

  if(dynamic){
    SIGi <- as.matrix(solve(Matrix::nearPD(SIG)$mat))  ####hessian
    mu4 <- as.vector(t(resid^2)%*%(resid^2)/(n*t)) #the 4th moment of residuals
    OMG<-f_OMG(sige,Gn,n,kz,mu4)
    varcov<-SIGi+SIGi%*%OMG%*%SIGi
    tmpplus<-diag(abs(varcov))
  }else {
    SIGi <- as.matrix(solve(Matrix::nearPD(SIG)$mat))
    varcov<-SIGi
    tmpplus <- diag(SIGi)
  }

  theta <-c(bhat,rho,sige)
  if(dynamic){  std <- sqrt(tmpplus)/sqrt(n*t)  }else {
    std <- sqrt(tmpplus)} ##standard errors
  tmps <- theta/std  ##t-statistic
  pval<-2*pnorm(abs(tmps),lower.tail=FALSE)

  ####Bias correction for non-dynamic model with twoway effect
  if(LYtrans==TRUE & effect == "twoways" & !dynamic) {
    bias01<-(-SIGi%*%c(rep(0,kz),1/(1-rho),1/(2*sige))/n)
    thetat<-theta-bias01
    if(!demn){
      bias02<-matrix(0,nrow=(kz+2),ncol=(kz+2))
      bias02[1:(kz+1),1:(kz+1)]<-diag(kz+1)
      bias02[kz+2,kz+2]<-t/(t-1)
      thetat<-bias02%*%thetat
    }
    names(thetat)<-names(theta)
    bhat<-thetat[1:kz]
    rho<-thetat[kz+1]
    sige<-thetat[kz+2]

    if(demn){ likl <-f2_sar2(rho,bhat,y0,Z0,detval,n0,t0,
                             sige,sind0,tind0,Wy0)}

    fsig2<-f_SIG(rho,bhat,sige,W,Z,n,t,kz)
    SIG<-fsig2$SIG; Gn<-fsig2$Gn #; Sni<-fsig2$Sni
    SIGi<-as.matrix(solve(SIG)); tmpplus<-diag(SIGi)
    varcov<-SIGi
    std <- sqrt(tmpplus)   ##standard errors
    tmps <- thetat/std ###tmps <- thetan/std  ##t-statistic
    pval<-2*pnorm(abs(tmps),lower.tail=FALSE) ###p-values
  }

  ####bias correction for dynamic panel
  if((dynamic & LYtrans & effect %in% c("individual","twoways")) ||
     (dynamic & DIRtrans & effect %in% c("twoways"))){
    iNm<-diag(n) ####iNm <- Matrix::Matrix(diag(n), sparse = TRUE)
    if(tlaginfo$tl & tlaginfo$stl){
      An<-Sni%*%(bhat[1]*iNm+bhat[2]*W)
    } else if(!tlaginfo$tl & tlaginfo$stl){
      An<-Sni%*%(bhat[1]*W)
    } else if(tlaginfo$tl & !tlaginfo$stl){
      An<-Sni%*%(bhat[1]*iNm)
    } else{ stop("Error dynamic structure!")}

    if(dynamic & LYtrans & effect %in% c("individual","twoways")){
      Rn<-Re(eigen(An)$vectors);   Dn<-Re(eigen(An)$values)
      Jn<-matrix(0,nrow = n,ncol = n)
      mm<-0; for(i in 1:n){ if(Dn[i]>1-1/n){Jn[i,i]<-Dn[i]; mm<-mm+1}}
      if(mm>=1){ Bn<-An-Rn%*%Jn%*%solve(Rn) }else { Bn<-An }
    }else if(dynamic & DIRtrans & effect %in% c("twoways")){
      Bn<-An
    }

    bias1s<-rep(0,kz+2);  bias1u<-rep(0,kz+2)
    if(tlaginfo$tl & tlaginfo$stl){
      bias1s[1] <- (1/n)*sum(diag(solve(iNm-Bn)%*%Sni))
      bias1s[2] <- (1/n)*sum(diag(W%*%solve(iNm-Bn)%*%Sni))
      bias1s[kz+1] <- (1/n)*sum(diag(Gn%*%solve(iNm-Bn)%*%Sni))*bhat[1]+
        (1/n)*sum(diag(Gn%*%W%*%solve(iNm-Bn)%*%Sni))*bhat[2]+
        (1/n)*sum(diag(Gn))
    } else if(!tlaginfo$tl & tlaginfo$stl){
      bias1s[1] <- (1/n)*sum(diag(solve(iNm-Bn)%*%Sni))
      bias1s[kz+1] <- (1/n)*sum(diag(Gn%*%solve(iNm-Bn)%*%Sni))*bhat[1]+
        (1/n)*sum(diag(Gn))
    } else if(tlaginfo$tl & !tlaginfo$stl){
      bias1s[1] <- (1/n)*sum(diag(solve(iNm-Bn)%*%Sni))
      bias1s[kz+1] <- (1/n)*sum(diag(Gn%*%solve(iNm-Bn)%*%Sni))*bhat[1]+
        (1/n)*sum(diag(Gn))
    }
    bias1s[kz+2]<-1/(2*sige)

    if(dynamic & LYtrans & effect %in% c("individual","twoways")){
      bias1utemp<-t/(2*(1-rho))
      if(tlaginfo$tl & tlaginfo$stl){
        bias1u[1]<-bias1utemp; bias1u[2]<-bias1utemp; bias1u[kz+1]<-bias1utemp
      }else if(!tlaginfo$tl & tlaginfo$stl){
        bias1u[1]<-bias1utemp; bias1u[kz+1]<-bias1utemp
      }else if(tlaginfo$tl & tlaginfo$stl){
        bias1u[1]<-bias1utemp; bias1u[kz+1]<-bias1utemp
      }
      if(mm>=1){  bias1<-bias1s+bias1u*(mm/n) }else{  bias1<-bias1s }
    } else{ bias1<-bias1s }

    bias<-(-SIGi%*%bias1/t)
    theta1<- theta-bias

    if(dynamic & DIRtrans & effect %in% c("twoways")){
      bias2<-matrix(0,nrow=(kz+2),ncol=1)
      bias2[kz+1,1]<-1/(1-rho)
      bias2[kz+2,1]<-1/(2*sige)
      bias_2<-(-SIGi%*%bias2/n)
      theta1<-theta1-bias_2
    }

    bhattemp<-theta1[1:kz]
    rhotemp<-theta1[kz+1]
    sigetemp<-theta1[kz+2]

    likl1<-f2_sar_dyn(rhotemp,bhattemp,y0,Z0,detval,
                      n0,t0,sigetemp,sind0,tind0,Wy0)
    ###f2_sar(rhotemp,bhattemp,y,Z,Wy,detval,n,t,sigetemp)

    fsig3<-f_SIG(rhotemp,bhattemp,sigetemp,W,Z,n,t,kz)
    SIGtemp<-fsig3$SIG; Gntemp<-fsig3$Gn
    SIGtemp <- SIGtemp/(n*t)
    SIGitemp <- solve(SIGtemp)

    yhat1 <- Matrix::solve(Matrix::Matrix(diag(n*t)- rhotemp*Wnt,
                                          sparse = TRUE))%*%(Z%*%bhattemp)
    yhat1 <- c(as.matrix(yhat1))
    resid1 <- y-yhat1

    mu4temp <- as.vector(t(resid1^2)%*%(resid1^2)/(n*t))
    OMGtemp<-f_OMG(sigetemp,Gntemp,n,kz,mu4temp)

    tmpplus1<-diag(abs(SIGitemp+SIGitemp%*%OMGtemp%*%SIGitemp))
    varcov<-(SIGitemp+SIGitemp%*%OMGtemp%*%SIGitemp)/(n*t)

    std1<-sqrt(tmpplus1)/sqrt(n*t)
    tmps1 <- theta1/std1
    pval1<-2*pnorm(abs(tmps1),lower.tail=FALSE)

    residr1 <- as.vector(y-rhotemp*Wy-Z%*%bhattemp)
    #ymean <- y0-mean(y0)
    ymean <- y-mean(y)
    rsqr2 <- crossprod(ymean)
    #rsqr1 <- crossprod(residr)
    rsqr1 <- as.vector(residr%*%residr1)
    rsqr <- 1-c(rsqr1/rsqr2)
    residuals<-resid1

    res1 <- as.vector(ymean)
    res2 <- as.vector(yhat1-mean(y))

  } else if(!demn & !dynamic){
    reseff<-feffects(rho,beta=bhat,as.numeric(sige),W0,y,X=Z,n0,t0,y0,X0=Z0,
                     mny,mnx,mty,mtx,effect,tind,sind,Wy0)
    ymean<-y0-mean(y0)
    rsqr2 <- crossprod(ymean)
    rsqr1 <- crossprod(reseff$res.e)
    rsqr <- 1-as.vector(rsqr1/rsqr2)
    residuals<-reseff$res.e

    res1 <- y-mean(y)
    res2 <- yhat-mean(y)
  } else {
    residr <- as.vector(y - rho*Wy - Z%*%bhat)
    ymean <- y0-mean(y0)
    rsqr2 <- crossprod(ymean)
    rsqr1 <- crossprod(residr)
    rsqr <- 1-c(rsqr1/rsqr2)
    residuals<-residr

    res1 <- y-mean(y)
    res2 <- yhat-mean(y)
  }

  rsq1 <- as.vector(as.vector(res1)%*%as.vector(res2))
  rsq2 <- as.vector(crossprod(res1))
  rsq3 <- as.vector(crossprod(res2))
  adjrsqr <- rsq1^2/(rsq2*rsq3)

  results<-list()
  results$coefficients<- c(bhat)
  results$rho <- rho
  results$rho.tst<-tmps[kz+1]
  results$rho.se<-std[kz+1]
  results$rho.pval<-pval[kz+1]
  results$tstat <- tmps[1:kz]
  results$std<-std[1:kz]
  results$pval<-pval[1:kz]
  results$sige<-sige
  results$likl <- likl
  if((dynamic & LYtrans & effect %in% c("individual","twoways")) ||
     (dynamic & DIRtrans & effect %in% c("twoways"))){
    results$coefficients1 <- theta1[1:kz]
    names(results$coefficients1)<-names(results$coefficients)
    results$rho1 <- theta1[kz+1]
    names(results$rho1)<-"rho1"
    results$rho.tst1<-tmps1[kz+1]
    results$rho.se1<-std1[kz+1]
    results$rho.pval1<-pval1[kz+1]
    results$tstat1 <- tmps1[1:kz]
    results$std1<-std1[1:kz]
    results$pval1<-pval1[1:kz]
    results$sige1<-theta1[kz+2]
    results$likl1<-likl1
  }
  results$rsqr<-rsqr
  results$adjrsqr<-adjrsqr
  results$varcov<-varcov
  results$effect<-effect
  results$model<-model
  results$call<-cl
  results$dynamic<-dynamic
  results$LeeYu<-LYtrans
  results$demn<-demn
  results$DirectT<-DIRtrans
  results$tlaginfo<-tlaginfo
  results$residuals<-residuals
  results$W<-W
  results$W0<-W0
  results$interval<-c(rmin,rmax)
  if(!demn & !DIRtrans & !dynamic){
    results$detval<-detval
    results$int.tab<-reseff$int.tab
    if(effect %in% c("time","twoways")){results$tfe.tab<-reseff$tfe.tab}
    if(effect %in% c("individual","twoways")){results$sfe.tab<-reseff$sfe.tab}
  }
  if(dynamic & tlaginfo$tl & tlaginfo$stl){
    if((dynamic & LYtrans & effect %in% c("individual","twoways")) ||
       (dynamic & DIRtrans & effect %in% c("twoways"))){
      res<-parWald(theta1,varcov)
    }else{res<-parWald(as.matrix(theta,ncol=1),varcov) }
    results$Waldt <- c(res$Waldt)
    results$pWald <- c(res$F1)
  }
  class(results) <- "SDPDm"
  return(results)
}




parWald<-function(theta1,varcov){
  npar<-length(theta1)
  tau <- theta1[1,1]
  eta <- theta1[2,1]
  rho <- theta1[npar-1,1]
  R<- tau+rho+eta
  Rafg<-rep(0,npar)
  Rafg[1]<-1;  Rafg[2]<-1;  Rafg[npar-1]<-1
  Waldt<-R*(solve(t(Rafg)%*%varcov%*%Rafg))*R
  F1<-1-pchisq(Waldt,1)
  result<-list(Waldt=Waldt,F1=F1)
  return(result)
}





