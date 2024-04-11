#' @name print.summary.SDPDm
#'
#' @title Print of summary for class SDPDm
#'
#' @description Method for printing the summary the results of objects of class "SDPDm"
#'
#' @method print summary.SDPDm
#'
#' @param x summary object of class "SDPDm"
#' @param ... additional arguments to be passed
#'
#' @return No return value
#'
#' @seealso
#' \code{SDPDm}
#'
#' @author Rozeta Simonovska
#' 
#' @export


print.summary.SDPDm <- function(x,...) {
  
    if(x$dynamic){
      cat(paste0(x$model," dynamic panel model with ",
                 x$effect, " fixed effects\n"))
    }else{
      cat(paste0(x$model," panel model with ",
                 x$effect, " fixed effects\n"))
    }
    cat("\nCall:\n")
    print(x$call)
    cat("\nSpatial autoregressive coefficient:\n")
    
    if((x$dynamic & x$LeeYu & x$effect %in%
        c("individual","twoways")) ||
       (x$dynamic & x$DirectT & x$effect %in% c("twoways"))){
      
      rtab <- cbind(x$rho1,x$rho.se1,x$rho.tst1,x$rho.pval1)
      rownames(rtab)<-"rho"
      colnames(rtab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
    }else{
      rtab <- cbind(x$rho,x$rho.se,x$rho.tst,x$rho.pval)
      rownames(rtab)<-"rho"
      colnames(rtab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
    }
    printCoefmat(rtab,  signif.legend=FALSE)
    cat("\nCoefficients:\n")
    if((x$dynamic & x$LeeYu & x$effect %in%
        c("individual","twoways")) ||
       (x$dynamic & x$DirectT & x$effect %in% c("twoways"))){
      ctab<-cbind(x$coefficients1,x$std1,x$tstat1,x$pval1)
      colnames(ctab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
      rownames(ctab) <- names(x$coefficients)
    }else{
      ctab<-cbind(x$coefficients,x$std,x$tstat,x$pval)
      colnames(ctab)<-c("Estimate","Std. Error","t-value","Pr(>|t|)")
      rownames(ctab) <- names(x$coefficients)
    }
    printCoefmat(ctab)
    cat("\n")
  
  invisible(x)
}
