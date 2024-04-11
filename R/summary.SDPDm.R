#' @name summary.SDPDm
#'
#' @title Summary for class SDPDm
#'
#' @description Method for summarizing the results of objects of class "SDPDm"
#'
#' @method summary SDPDm
#'
#' @param object object of class "SDPDm"
#' @param ... additional arguments to be passed
#'
#' @return Summary of SDPDm
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#' 
#' @export


summary.SDPDm <- function(object,...) {
  if(inherits(object,"SDPDm")){
    
    if((object$dynamic & object$LeeYu & object$effect %in%
        c("individual","twoways")) ||
       (object$dynamic & object$DirectT & object$effect %in% c("twoways"))){

      rtab <- cbind(object$rho1,object$rho.se1,object$rho.tst1,object$rho.pval1)
      rownames(rtab)<-"rho"
      colnames(rtab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
    }else{
      rtab <- cbind(object$rho,object$rho.se,object$rho.tst,object$rho.pval)
      rownames(rtab)<-"rho"
      colnames(rtab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
    }

    if((object$dynamic & object$LeeYu & object$effect %in%
        c("individual","twoways")) ||
       (object$dynamic & object$DirectT & object$effect %in% c("twoways"))){
      ctab<-cbind(object$coefficients1,object$std1,object$tstat1,object$pval1)
      colnames(ctab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
      rownames(ctab) <- names(object$coefficients)
    }else{
      ctab<-cbind(object$coefficients,object$std,object$tstat,object$pval)
      colnames(ctab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
      rownames(ctab) <- names(object$coefficients)
    }
    
    object$spatialcoefTab<-rtab
    object$coefTab<-ctab

    class(object) <- c("summary.SDPDm","SDPDm")
    object
  }

}
