#' @name impactsSDPDm
#'
#' @title Impacts for 'SDPDm' objects
#'
#' @description Direct and indirect effects estimates
#'
#' @param res an object of class 'SDPDm'
#' @param NSIM number of simulations to be performed, default = 200
#' @param sd starting seed, default = 12345
#'
#' @return An object of class 'impactsSDPDm'
#'
#' @details
#' For dynamic panel data model:
#' \deqn{y_{t} = \tau y_{t-1} + \rho W y_{t} + \eta W y_{t-1} + X_{t} \beta + W X_{t} \theta + \alpha + \mu + u_{t}}
#' Short term effects for k\emph{th} explanatory variable:
#' \deqn{(I - \rho W)^{-1}(\beta_{k} I_{n} + \theta_{k} W)}
#' Long term effects for k\emph{th} explanatory variable:
#' \deqn{((1-\tau)I_{n} - (\rho+\eta)W)^{-1}(\beta_{k} I_{n} + \theta_{k} W)}
#'
#' @author Rozeta Simonovska
#'
#' @seealso \code{\link{SDPDm}}
#'
#' @export

impactsSDPDm<-function(res, NSIM = 200, sd = 12345){
  if(!inherits(res,"SDPDm")){ stop("Wrong class object!")}

  npar <- length(res$coefficients)
  Wmat<-res$W0
  N<-nrow(Wmat)

  if(res$dynamic){
    if(res$tlaginfo$tl & res$tlaginfo$stl){
      px <- npar-2
    } else if(!res$tlaginfo$tl & res$tlaginfo$stl){
      px <- npar-1
    }else if(res$tlaginfo$tl & !res$tlaginfo$stl){
      px <- npar-1
    }else stop("ERROR in dynamic values!")
    if(res$model=="sdm"){px<-px/2}
    if((res$LeeYu & res$effect %in% c("individual","twoways")) ||
       (res$DirectT & res$effect %in% c("twoways"))){
      varb<-c(res$coefficients1,res$rho1,res$sige1)
    }else {varb<-c(res$coefficients,res$rho,res$sige)}
  }else{
    px <- npar
    if(res$model=="sdm"){px<-px/2}
    varb<-c(res$coefficients,res$rho,res$sige)
  }

  simdirst <- matrix(0,px,NSIM)
  simindst <- matrix(0,px,NSIM)
  simtotst <- matrix(0,px,NSIM)
  simdirlt <- matrix(0,px,NSIM)
  simindlt <- matrix(0,px,NSIM)
  simtotlt <- matrix(0,px,NSIM)

  for(sim in 1:NSIM){
    set.seed(sd+10*sim)
    parms <- t(chol(res$varcov))%*%rnorm(npar+2) + varb
    rhosim <- parms[npar+1,1] ##rho coef
    if(res$dynamic){
      if(res$tlaginfo$tl & res$tlaginfo$stl){
        betasim <- parms[3:npar,1]
        tausim <- parms[1,1]
        etasim <- parms[2,1]
      } else if(!res$tlaginfo$tl & res$tlaginfo$stl){
        betasim <- parms[2:npar,1]
        etasim <- parms[1,1]
      }else if(res$tlaginfo$tl & !res$tlaginfo$stl){
        betasim <- parms[2:npar,1]
        tausim <- parms[1,1]
      }else stop("ERROR in dynamic values!")
    } else {
      betasim <- parms[1:npar,1]
    }
    SS<-solve((diag(N)-rhosim*Wmat))

    for(p in 1:px){
      C <- matrix(0,N,N)
      if(res$model=="sdm") {  C<-betasim[px+p]*Wmat }
      diag(C)<-rep(betasim[p],N)

      SC <- SS%*%C
      simdirst[p,sim] <- sum(diag(SC))/N # average direct effect
      simtotst[p,sim] <- sum(SC)/N
      simindst[p,sim] <- simtotst[p,sim] - simdirst[p,sim]

      if(res$dynamic){
        if(res$tlaginfo$tl & res$tlaginfo$stl){
          SC <- solve((1-tausim)*diag(N)-(rhosim+etasim)*Wmat)%*%C
        }else if(!res$tlaginfo$tl & res$tlaginfo$stl){
          SC <- solve(diag(N)-(rhosim+etasim)*Wmat)%*%C
        }else if(res$tlaginfo$tl & !res$tlaginfo$stl){
          SC <- solve((1-tausim)*diag(N)-(rhosim)*Wmat)%*%C
        }
        simdirlt[p,sim] <- sum(diag(SC))/N   # average direct effect
        simtotlt[p,sim] <- sum(SC)/N
        simindlt[p,sim] <- simtotlt[p,sim] - simdirlt[p,sim]
      }
    }
  }

  ###Short-term direct effects
  DirSt<-rowMeans(simdirst)
  sdDSt<-apply(simdirst,1,sd)
  tstDSt<-DirSt/sdDSt
  pvalDSt<-2*pnorm(abs(tstDSt),lower.tail=FALSE)
  ###Short-term indirect effects
  IdirSt<-rowMeans(simindst)
  sdISt<-apply(simindst,1,sd)
  tstISt<-IdirSt/sdISt
  pvalISt<-2*pnorm(abs(tstISt),lower.tail=FALSE)
  ###Short-term total effects
  TotSt<-rowMeans(simtotst)
  sdTSt<-apply(simtotst,1,sd)
  tstTSt<-TotSt/sdTSt
  pvalTSt<-2*pnorm(abs(tstTSt),lower.tail=FALSE)

  if(res$dynamic){
    if(res$tlaginfo$tl & res$tlaginfo$stl){
      vnam<-c(names(res$coefficients)[-1])
      vnam<-vnam[-1]
    }else if((!res$tlaginfo$tl & res$tlaginfo$stl) ||
             (res$tlaginfo$tl & !res$tlaginfo$stl)){
      vnam <- c(names(res$coefficients)[-1])
    }

    ###long-term direct effects
    DirLt<-rowMeans(simdirlt)
    sdDLt<-apply(simdirlt,1,sd)
    tstDLt<-DirLt/sdDLt
    pvalDLt<-2*pnorm(abs(tstDLt),lower.tail=FALSE)
    ###long-term indirect effects
    IdirLt<-rowMeans(simindlt)
    sdILt<-apply(simindlt,1,sd)
    tstILt<-IdirLt/sdILt
    pvalILt<-2*pnorm(abs(tstILt),lower.tail=FALSE)
    ###long-term total effects
    TotLt<-rowMeans(simtotlt)
    sdTLt<-apply(simtotlt,1,sd)
    tstTLt<-TotLt/sdTLt
    pvalTLt<-2*pnorm(abs(tstTLt),lower.tail=FALSE)

  }else{   vnam<- c(names(res$coefficients)) }
  if(res$model=="sdm"){ vnam<-vnam[1:px]  }

  result<-list()
  ###Direct short term
  DIRECTst.tab <- cbind(DirSt,sdDSt,tstDSt,pvalDSt)
  colnames(DIRECTst.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  rownames(DIRECTst.tab)<-vnam
  result$DIRECTst.tab<-DIRECTst.tab
  ###indirect short term
  INDIRECTst.tab <- cbind(IdirSt,sdISt,tstISt,pvalISt)
  colnames(INDIRECTst.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  rownames(INDIRECTst.tab)<-vnam
  result$INDIRECTst.tab<-INDIRECTst.tab
  ###total short term
  TOTALst.tab <- cbind(TotSt,sdTSt,tstTSt,pvalTSt)
  colnames(TOTALst.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  rownames(TOTALst.tab)<-vnam
  result$TOTALst.tab<-TOTALst.tab

  if(res$dynamic){
    ###Direct long term
    DIRECTlt.tab <- cbind(DirLt,sdDLt,tstDLt,pvalDLt)
    colnames(DIRECTlt.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    rownames(DIRECTlt.tab)<-vnam
    result$DIRECTlt.tab<-DIRECTlt.tab
    ###Indirect long term
    INDIRECTlt.tab <- cbind(IdirLt,sdILt,tstILt,pvalILt)
    colnames(INDIRECTlt.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    rownames(INDIRECTlt.tab)<-vnam
    result$INDIRECTlt.tab<-INDIRECTlt.tab
    ###Total long term
    TOTALlt.tab <- cbind(TotLt,sdTLt,tstTLt,pvalTLt)
    colnames(TOTALlt.tab) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
    rownames(TOTALlt.tab)<-vnam
    result$TOTALlt.tab<-TOTALlt.tab
  }
  class(result) <- c("impactsSDPDm")
  return(result)
}
